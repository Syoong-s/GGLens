#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include "esd_plot.hpp"
#include "esd_plot_source.inc"

namespace bcg_esd {
namespace {
using PyHandle = std::unique_ptr<PyObject, decltype(&Py_DecRef)>;

// ==========================================
// Function: checked
// Method: Adopt a Python reference or surface the active Python exception.
// ==========================================
PyHandle checked(PyObject *object) {
    if (!object) {
        PyErr_Print();
        throw std::runtime_error("In-process plotting failed; see Python traceback");
    }
    return PyHandle(object, &Py_DecRef);
}

// ==========================================
// Function: fd_plot_data
// Method: Marshal both Cartesian FD fits and their bin measurements with no text serialization.
// ==========================================
PyHandle fd_plot_data(const Config &cfg, const SelfCalibrationDiagnostics &diagnostics) {
    if (!cfg.self_calibrate) {
        Py_INCREF(Py_None);
        return checked(Py_None);
    }
    auto components = checked(PyTuple_New(2));
    for (int i = 0; i < 2; ++i) {
        const auto &points = diagnostics.points[i];
        const auto &fit = diagnostics.fits[i];
        auto rows = checked(PyList_New(static_cast<Py_ssize_t>(points.size())));
        for (std::size_t j = 0; j < points.size(); ++j) {
            const auto &p = points[j];
            auto row = checked(Py_BuildValue("(Ldddi)", static_cast<long long>(p.count),
                p.field, p.value, p.sigma, p.status));
            PyList_SET_ITEM(rows.get(), static_cast<Py_ssize_t>(j), row.release());
        }
        auto component = checked(Py_BuildValue("(ddddiO)", fit.m, fit.c, fit.sigma_m,
            fit.sigma_c, fit.count, rows.get()));
        PyTuple_SET_ITEM(components.get(), i, component.release());
    }
    return checked(Py_BuildValue("(sLdO)", cfg.fd_plot.c_str(),
        static_cast<long long>(diagnostics.selected), cfg.fd_limit(), components.get()));
}

// ==========================================
// Function: call_renderer
// Method: Load the bundled module and send the rank-zero measurement arrays to it.
// ==========================================
void call_renderer(const std::vector<Result> &results, const Config &cfg,
                   const SelfCalibrationDiagnostics &diagnostics) {
    auto module = checked(PyModule_New("esd_plot"));
    PyObject *globals = PyModule_GetDict(module.get());
    auto loaded = checked(PyRun_String(plot_source, Py_file_input, globals, globals));
    auto render = checked(PyObject_GetAttrString(module.get(), "render_plots"));
    auto rows = checked(PyList_New(static_cast<Py_ssize_t>(results.size())));
    for (std::size_t i = 0; i < results.size(); ++i) {
        const auto &r = results[i];
        auto row = checked(Py_BuildValue("(ddddddd)", r.effective, r.measured, r.sigma, r.cross,
            r.cross_sigma, r.theory, r.ratio));
        PyList_SET_ITEM(rows.get(), static_cast<Py_ssize_t>(i), row.release());
    }
    auto limits = cfg.plot_ratio_min != 0 || cfg.plot_ratio_max != 0
        ? checked(Py_BuildValue("(dd)", cfg.plot_ratio_min, cfg.plot_ratio_max))
        : checked(Py_BuildValue(""));
    auto fd = fd_plot_data(cfg, diagnostics);
    auto result = checked(PyObject_CallFunction(render.get(), "OssddOO", rows.get(),
        cfg.esd_plot.c_str(), cfg.plot_title.c_str(), cfg.plot_x_spread, cfg.plot_linear_width,
        limits.get(), fd.get()));
}
} // namespace

// ==========================================
// Function: render_plots
// Method: Initialize the linked Python environment locally for one end-of-run plotting call.
// ==========================================
void render_plots(const std::vector<Result> &results, const Config &cfg,
                  const SelfCalibrationDiagnostics &diagnostics) {
    PyConfig python;
    PyConfig_InitIsolatedConfig(&python);
    python.install_signal_handlers = 0;
    python.write_bytecode = 0;
    auto status = PyConfig_SetBytesString(&python, &python.program_name, FQ_PLOT_PYTHON);
    if (!PyStatus_Exception(status)) status = Py_InitializeFromConfig(&python);
    PyConfig_Clear(&python);
    if (PyStatus_Exception(status))
        throw std::runtime_error(std::string("Cannot initialize embedded Python: ") +
                                 (status.err_msg ? status.err_msg : "unknown error"));
    try {
        call_renderer(results, cfg, diagnostics);
    } catch (...) {
        Py_FinalizeEx();
        throw;
    }
    if (Py_FinalizeEx() < 0) throw std::runtime_error("Embedded Python finalization failed");
}
} // namespace bcg_esd
