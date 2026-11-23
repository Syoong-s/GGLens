#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <cstring>

#include "esd_fit_backend.hpp"
#include "esd_nfw_source.inc"
#include "esd_two_halo_source.inc"
#include "esd_plot_source.inc"
#include "esd_profiles_source.inc"
#include "esd_fit_source.inc"

namespace bcg_esd {
namespace {
using Object = std::unique_ptr<PyObject, decltype(&Py_DecRef)>;

// ==========================================
// Function: checked
// Method: Adopt a Python result and print the traceback on failure.
// ==========================================
Object checked(PyObject *raw) {
    if (!raw) {
        PyErr_Print();
        throw std::runtime_error("Embedded ESD fit failed; see Python traceback");
    }
    return Object(raw, &Py_DecRef);
}

// ==========================================
// Function: load_module
// Method: Import a compiled-in Python module by name for normal intra-library imports.
// ==========================================
Object load_module(const char *name, const char *source) {
    auto module = checked(PyModule_New(name));
    if (PyDict_SetItemString(PyImport_GetModuleDict(), name, module.get()) != 0)
        checked(nullptr);
    auto code = checked(Py_CompileString(source, name, Py_file_input));
    auto *globals = PyModule_GetDict(module.get());
    checked(PyEval_EvalCode(code.get(), globals, globals));
    return module;
}

// ==========================================
// Function: doubles
// Method: Unpack a contiguous native float64 buffer into MPI storage.
// ==========================================
std::vector<double> doubles(PyObject *bytes) {
    char *buffer = nullptr;
    Py_ssize_t size = 0;
    if (PyBytes_AsStringAndSize(bytes, &buffer, &size) != 0) checked(nullptr);
    if (size % sizeof(double)) throw std::runtime_error("Misaligned fitter float64 buffer");
    std::vector<double> values(size / sizeof(double));
    std::memcpy(values.data(), buffer, size);
    return values;
}

// ==========================================
// Function: packed
// Method: Give NumPy a byte buffer without building a separate object per sample.
// ==========================================
Object packed(const std::vector<double> &values) {
    return checked(PyBytes_FromStringAndSize(reinterpret_cast<const char *>(values.data()),
                                            values.size()*sizeof(double)));
}

// ==========================================
// Function: as_list
// Method: Transfer the small radial and covariance arrays as numeric sequences.
// ==========================================
Object as_list(const std::vector<double> &values) {
    auto list = checked(PyList_New(values.size()));
    for (std::size_t i = 0; i < values.size(); ++i)
        PyList_SET_ITEM(list.get(), i, checked(PyFloat_FromDouble(values[i])).release());
    return list;
}

// ==========================================
// Class: Interpreter
// Method: Use the linked CPython installation and leave MPI signal handling in C++.
// ==========================================
class Interpreter {
public:
    Interpreter() {
        PyConfig config;
        PyConfig_InitIsolatedConfig(&config);
        config.install_signal_handlers = 0;
        config.write_bytecode = 0;
        auto status = PyConfig_SetBytesString(&config, &config.program_name, FQ_ESD_PYTHON);
        if (!PyStatus_Exception(status)) status = Py_InitializeFromConfig(&config);
        const bool failed = PyStatus_Exception(status);
        const std::string message = failed && status.err_msg ? status.err_msg : "unknown error";
        PyConfig_Clear(&config);
        if (failed) throw std::runtime_error("Cannot initialize embedded Python: " + message);
    }

    ~Interpreter() { Py_FinalizeEx(); }
};
} // namespace

struct FitBackend::Runtime {
    Interpreter python;
    Object module{nullptr, &Py_DecRef};
    Object problem{nullptr, &Py_DecRef};
};

// ==========================================
// Function: FitBackend
// Method: Load the compiled-in model and use an optional isolated dependency path.
// ==========================================
FitBackend::FitBackend(const Config &cfg) : runtime_(std::make_unique<Runtime>()) {
    if (!cfg.fit_python_path.empty()) {
        auto path = checked(PyUnicode_DecodeFSDefault(cfg.fit_python_path.c_str()));
        if (PyList_Insert(PySys_GetObject("path"), 0, path.get()) != 0) checked(nullptr);
    }
    load_module("esd_nfw", esd_nfw_source);
    load_module("esd_two_halo", esd_two_halo_source);
    load_module("esd_profiles", esd_profiles_source);
    load_module("esd_plot", esd_plot_source);
    runtime_->module = load_module("esd_fit", esd_fit_source);
}

FitBackend::~FitBackend() = default;

// ==========================================
// Function: correlation
// Method: Generate the nonlinear matter correlation once for the whole stack.
// ==========================================
std::vector<double> FitBackend::correlation(const Config &cfg, const LensStack &stack) {
    auto corr = checked(PyObject_CallMethod(runtime_->module.get(), "nonlinear_correlation", "ddddddd",
                                           cfg.omega_m, cfg.h, cfg.omega_b_h2, cfg.spectral_index,
                                           cfg.sigma8, cfg.neutrino_mass, stack.redshift));
    return doubles(corr.get());
}

// ==========================================
// Function: setup
// Method: Provide the measurement and covariance arrays directly to the likelihood.
// ==========================================
void FitBackend::setup(const FitInput &input, const std::vector<double> &corr,
                       const Config &cfg, const LensStack &stack) {
    const auto radius = as_list(input.radius);
    const auto measured = as_list(input.measured);
    const auto covariance = as_list(input.covariance);
    const auto data = packed(corr);
    runtime_->problem = checked(PyObject_CallMethod(runtime_->module.get(), "FitProblem", "OOOddO",
                               radius.get(), measured.get(), covariance.get(),
                               stack.stellar_mass, cfg.omega_m, data.get()));
}

// ==========================================
// Function: sample
// Method: Run one independent ensemble on the current rank, retaining post-burn draws.
// ==========================================
FitChain FitBackend::sample(const Config &cfg, int ensemble) {
    auto result = checked(PyObject_CallMethod(runtime_->module.get(), "sample_ensemble", "Oiiiii",
                                              runtime_->problem.get(), cfg.fit_walkers, cfg.fit_steps,
                                              cfg.fit_burn, cfg.fit_seed, ensemble));
    return {doubles(PyTuple_GetItem(result.get(), 0)), doubles(PyTuple_GetItem(result.get(), 1))};
}

// ==========================================
// Function: finish
// Method: Pass product paths and metadata for final posterior export and rendering.
// ==========================================
void FitBackend::finish(const FitChain &chain, const Config &cfg, const LensStack &stack,
                        const std::vector<int> &fit_bins) {
    const auto samples = packed(chain.samples);
    const auto stats = packed(chain.diagnostics);
    const auto shape = checked(Py_BuildValue("(iii)", cfg.fit_ensembles,
                                            cfg.fit_steps - cfg.fit_burn, cfg.fit_walkers));
    const auto paths = checked(Py_BuildValue("{s:s,s:s,s:s,s:s,s:s}",
        "summary", cfg.summary_output.c_str(), "chain", cfg.chain_output.c_str(),
        "profiles", cfg.profile_output.c_str(), "draws", cfg.profile_draws_output.c_str(),
        "plot", cfg.fit_plot.c_str()));
    const auto bins = checked(PyList_New(fit_bins.size()));
    for (std::size_t i = 0; i < fit_bins.size(); ++i)
        PyList_SET_ITEM(bins.get(), i, checked(PyLong_FromLong(fit_bins[i])).release());
    const auto context = checked(Py_BuildValue("{s:s,s:d,s:n,s:i,s:i,s:i,s:i,s:s,s:O}",
        "center", bcg::center_name(cfg.center), "redshift", stack.redshift,
        "n_lenses", static_cast<Py_ssize_t>(stack.n_lenses), "regions", stack.regions,
        "burn_in", cfg.fit_burn, "profile_draws", cfg.fit_profile_draws, "seed", cfg.fit_seed,
        "title", cfg.plot_title.c_str(), "fit_bins", bins.get()));
    checked(PyObject_CallMethod(runtime_->module.get(), "finish_fit", "OOOOOO",
                               runtime_->problem.get(), samples.get(), stats.get(),
                               shape.get(), paths.get(), context.get()));
}
} // namespace bcg_esd
