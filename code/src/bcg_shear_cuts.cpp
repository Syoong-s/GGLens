#include "bcg_shear_cuts.hpp"

namespace bcg {
namespace {
using Field = ShearRow::Field;

// ==========================================
// Function: select_shear_row
// Method: Apply the active cuts and prepare the measurements and optional stellar-bar inputs.
// ==========================================
bool select_shear_row(const ShearRow &row, const Cuts &cuts, Source &source, BarDatum &bar) {
    if (!row.passes_column_cuts) return false;
    const auto &v = row.values;
    // Use the nearest integer CCD ID, as in the FD catalog filter.
    if (cuts.ccd_cut_enabled() &&
        std::binary_search(cuts.excluded_ccds.begin(), cuts.excluded_ccds.end(),
                           std::round(v[Field::ccd_num]))) return false;

    const double ra = v[Field::ra], dec = v[Field::dec], zp = v[Field::zp];
    const double gf1 = v[Field::gf1], gf2 = v[Field::gf2];

    // Sky coordinates and photometric selection.
    if (ra < 0 || ra >= 360 || dec < -90 || dec > 90) return false;
    if (cut_enabled(cuts.zplow) && zp <= cuts.zplow) return false;
    if (cut_enabled(cuts.zphigh) && zp >= cuts.zphigh) return false;

    // Exposure, PSF, and stamp quality.
    if (cut_enabled(cuts.chi2_thresh) && v[Field::chi2] > cuts.chi2_thresh) return false;
    if (cut_enabled(cuts.psf_chi2_mltp) && v[Field::poly_chi2] > cuts.psf_chi2_mltp) return false;
    if (cut_enabled(cuts.starcut) && v[Field::nstar] < cuts.starcut) return false;
    if (cut_enabled(cuts.flagcut) && v[Field::flag] * v[Field::sigma] > cuts.flagcut) return false;
    if (cut_enabled(cuts.imaxcut) && v[Field::imax] >= cuts.imaxcut) return false;
    if (cut_enabled(cuts.jmaxcut) && v[Field::jmax] >= cuts.jmaxcut) return false;
    if (cut_enabled(cuts.gf_lim) && (std::abs(gf1) > cuts.gf_lim || std::abs(gf2) > cuts.gf_lim)) return false;

    // SNR and size use half-light area. The size ratio also uses PSF FWHM.
    const bool need_snr = cuts.needs_snr(), need_psf = cuts.needs_psf();
    const double area = v[Field::half_light_area], psf = v[Field::psf_fwhm];
    if ((need_snr || need_psf) && area <= 0) return false;
    if (need_psf && psf <= 0) return false;
    if (cut_enabled(cuts.snrfcut) && v[Field::snr_f] < cuts.snrfcut) return false;
    const double snr = need_snr ? v[Field::half_light_flux] / std::sqrt(area) : 0.0;
    if (cut_enabled(cuts.snrlow) && snr < cuts.snrlow) return false;
    if (cut_enabled(cuts.snrhigh) && snr > cuts.snrhigh) return false;
    constexpr double pi = 3.14159265358979323846;
    const double rhalf = need_psf ? std::sqrt(area / pi) * 0.2628 * 2.0 : 0.0;
    if (cut_enabled(cuts.r_half_thresh) && rhalf <= cuts.r_half_thresh * psf) return false;

    source = Source({ra, dec, zp, v[Field::g1], v[Field::g2], v[Field::de],
                     v[Field::h1], v[Field::h2], gf1, gf2, v[Field::zp_err]});
    if (cuts.bar_enabled()) {
        bar = {v[Field::mag_stellar_locus], (rhalf - psf) / psf, snr, 0};
        if (cuts.exposure_bar_enabled()) {
            const double exposure = v[Field::expo_num];
            if (exposure < 1 || exposure > 10000000)
                throw std::runtime_error("Invalid exposure ID for enabled exposure bar");
            bar.exposure_id = static_cast<int>(std::llround(exposure));
            if (std::abs(exposure - bar.exposure_id) > 0.01)
                throw std::runtime_error("Non-integer exposure ID for enabled exposure bar");
        }
    }
    return true;
}

// ==========================================
// Function: filter_duplicate_groups
// Method: Reject over-limit adjacent sky-position groups using a bounded row buffer.
// ==========================================
void filter_duplicate_groups(const std::string &path, const ShearLayout &layout, std::size_t max_dup,
                             const std::function<void(const ShearRow &)> &visit) {
    constexpr double coordinate_tolerance = 1e-7;
    std::size_t count = 0;
    double first_ra = 0, first_dec = 0;
    std::vector<ShearRow> group;
    group.reserve(std::min<std::size_t>(max_dup, 64));
    const auto flush = [&] {
        if (count <= max_dup)
            for (const auto &row : group) visit(row);
        group.clear();
        count = 0;
    };
    for_each_shear_row(path, layout, [&](const ShearRow &row) {
        const double ra = row.values[Field::ra], dec = row.values[Field::dec];
        if (count && (std::abs(ra - first_ra) > coordinate_tolerance ||
                      std::abs(dec - first_dec) > coordinate_tolerance)) flush();
        if (count == 0) {
            first_ra = ra;
            first_dec = dec;
        }
        ++count;
        if (count <= max_dup) group.push_back(row);
    });
    flush();
}
} // namespace

// ==========================================
// Function: for_each_local_source
// Method: Select worker-owned rows using the duplicate and per-row cuts.
// ==========================================
void for_each_local_source(const std::vector<std::string> &paths, const Cuts &cuts,
                           const ShearLayout &layout, MPI_Comm comm,
                           const std::function<void(const Source &, const BarDatum &)> &visit) {
    int rank = 0, size = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    const auto select = [&](const ShearRow &row) {
        Source source;
        BarDatum bar;
        if (select_shear_row(row, cuts, source, bar)) visit(source, bar);
    };
    for (std::size_t i = 0; i < paths.size(); ++i) {
        if (size > 1 && (rank == 0 || static_cast<int>(i % (size - 1)) != rank - 1)) continue;
        if (cut_enabled(cuts.max_dup))
            filter_duplicate_groups(paths[i], layout, static_cast<std::size_t>(cuts.max_dup), select);
        else
            for_each_shear_row(paths[i], layout, select);
    }
}

// ==========================================
// Function: read_local_sources
// Method: Collect selected rows and their optional aligned stellar-bar metadata.
// ==========================================
LocalSources read_local_sources(const std::vector<std::string> &paths, const Cuts &cuts,
                                const ShearLayout &layout, MPI_Comm comm) {
    LocalSources result;
    const bool need_bar = cuts.bar_enabled();
    for_each_local_source(paths, cuts, layout, comm, [&](const Source &source, const BarDatum &bar) {
        result.sources.push_back(source);
        if (need_bar) result.bar_data.push_back(bar);
    });
    return result;
}
} // namespace bcg
