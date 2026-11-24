#include "bcg_archive.hpp"
#include "bcg_esd.hpp"
#include <iostream>
#include <numeric>

namespace bcg_esd {
namespace {
constexpr char magic[] = "FQESDMEA";
constexpr std::uint32_t version = 1;

std::string calibration_summary(const Measurement &measurement) {
    std::ostringstream out;
    if (measurement.calibration_sources >= 0)
        out << "Self-calibration all-R unique rows=" << measurement.calibration_sources << "; ";
    out << measurement.calibration.description();
    return out.str();
}
} // namespace

bool is_measurement_snapshot(const bcg::fs::path &path) {
    return bcg::archive::has_magic(path, magic);
}

void write_measurement(const bcg::fs::path &path, const Measurement &m) {
    bcg::archive::publish(path, [&](const auto &temporary) {
        bcg::archive::Writer out(temporary);
        const auto &[center, cosmology, sample, radial_bins] = m.metadata;
        const auto &cal = m.calibration;
        out.write(magic, sizeof(magic) - 1);
        out.put(version);
        out.put(static_cast<std::uint32_t>(center));
        out.put(static_cast<std::uint32_t>(radial_bins));
        out.put(static_cast<std::uint32_t>(sample.regions));
        out.put(sample.redshift);
        out.put(sample.stellar_mass);
        const double cosmology_values[] = {cosmology.omega_m, cosmology.h, cosmology.sigma8,
            cosmology.omega_b_h2, cosmology.spectral_index, cosmology.neutrino_mass};
        out.write(cosmology_values, std::size(cosmology_values));
        out.put(m.sources_after_cuts);
        out.put(m.calibration_sources);
        out.put(m.photometry_selected);
        out.put(static_cast<std::uint32_t>(m.requested_regions));
        out.put(static_cast<std::uint32_t>(m.fd_restored));
        out.put(static_cast<std::uint32_t>(m.fd_bins));
        out.put(m.fd_limit);
        const double calibration_values[] = {cal.m1, cal.m2, cal.c1, cal.c2};
        out.write(calibration_values, std::size(calibration_values));

        for (int groups : m.groups_per_region) out.put(static_cast<std::uint32_t>(groups));
        out.put(m.bins.front().low);
        for (const auto &bin : m.bins) out.put(bin.high);
        for (const auto &bin : m.bins) {
            out.put(bin.effective);
            out.put(bin.count);
            out.put(bin.measured);
            out.put(bin.cross);
        }
        for (const auto &bin : m.bins) out.write(bin.tangential_replicas.data(), sample.regions);
        for (const auto &bin : m.bins) out.write(bin.cross_replicas.data(), sample.regions);
        out.seal();
    });
}

Measurement read_measurement(const bcg::fs::path &path) {
    try {
        bcg::archive::Reader in(path);
        char prefix[sizeof(magic) - 1];
        in.read(prefix, sizeof(prefix));
        if (std::string(prefix, sizeof(prefix)) != magic || in.get<std::uint32_t>() != version)
            throw std::runtime_error("Unsupported measurement snapshot format");
        const auto center = in.get<std::uint32_t>();
        const auto n = in.get<std::uint32_t>();
        const auto regions = in.get<std::uint32_t>();
        if (center > static_cast<std::uint32_t>(bcg::Center::nwc) || n < 1 || n > max_radial_bins ||
            regions < 2 || regions >= std::numeric_limits<int>::max())
            throw std::runtime_error("Invalid measurement dimensions");

        Measurement m;
        auto &sample = m.metadata.sample;
        m.metadata.center = static_cast<bcg::Center>(center);
        m.metadata.radial_bins = static_cast<int>(n);
        sample.regions = static_cast<int>(regions);
        sample.redshift = in.get<double>();
        sample.stellar_mass = in.get<double>();
        double cosmology[6];
        in.read(cosmology, std::size(cosmology));
        m.metadata.cosmology = {cosmology[0], cosmology[1], cosmology[2], cosmology[3], cosmology[4], cosmology[5]};
        m.sources_after_cuts = in.get<std::int64_t>();
        m.calibration_sources = in.get<std::int64_t>();
        m.photometry_selected = in.get<std::int64_t>();
        const auto requested = in.get<std::uint32_t>();
        const auto restored = in.get<std::uint32_t>();
        const auto fd_bins = in.get<std::uint32_t>();
        m.fd_limit = in.get<double>();
        double calibration[4];
        in.read(calibration, std::size(calibration));
        m.calibration = {calibration[0], calibration[1], calibration[2], calibration[3]};
        if (requested < regions || requested > std::numeric_limits<int>::max() || restored > 1 ||
            fd_bins > 100 || m.sources_after_cuts <= 0 || m.calibration_sources < -1 ||
            m.calibration_sources > m.sources_after_cuts || m.photometry_selected < -1 ||
            (fd_bins == 0 ? m.calibration_sources != -1 || m.fd_limit != 0 :
             fd_bins < 2 || m.calibration_sources <= 0 || !(m.fd_limit > 0 && std::isfinite(m.fd_limit))))
            throw std::runtime_error("Invalid measurement statistics");
        m.requested_regions = static_cast<int>(requested);
        m.fd_restored = restored != 0;
        m.fd_bins = static_cast<int>(fd_bins);
        validate_cosmology(m.metadata.cosmology);
        m.calibration.validate();

        const std::uint64_t expected = 4ULL*regions + 8ULL*(2 + 5ULL*n + 2ULL*n*regions);
        if (in.remaining() != expected) throw std::runtime_error("Invalid measurement snapshot length");
        m.groups_per_region.resize(regions);
        for (auto &groups : m.groups_per_region) {
            const auto count = in.get<std::uint32_t>();
            if (!count || count > std::numeric_limits<int>::max())
                throw std::runtime_error("Invalid jackknife region group count");
            groups = static_cast<int>(count);
            sample.n_lenses += groups;
        }
        if (m.photometry_selected > sample.n_lenses ||
            (m.metadata.center != bcg::Center::bcg && m.photometry_selected != -1))
            throw std::runtime_error("Invalid photometry selection count");
        m.bins.resize(n);
        auto edge = in.get<double>();
        for (auto &bin : m.bins) {
            bin.low = edge;
            bin.high = edge = in.get<double>();
            if (!(bin.low > 0 && bin.high >= bin.low && std::isfinite(bin.high)))
                throw std::runtime_error("Invalid radial edges");
        }
        for (auto &bin : m.bins) {
            bin.effective = in.get<double>();
            bin.count = in.get<std::int64_t>();
            bin.measured = in.get<double>();
            bin.cross = in.get<double>();
            if (bin.count < 0) throw std::runtime_error("Invalid lens–source pair count");
        }
        for (auto &bin : m.bins) {
            bin.tangential_replicas.resize(regions);
            in.read(bin.tangential_replicas.data(), regions);
        }
        for (auto &bin : m.bins) {
            bin.cross_replicas.resize(regions);
            in.read(bin.cross_replicas.data(), regions);
        }
        m.id = in.verify();
        return m;
    } catch (const std::exception &error) {
        throw std::runtime_error("Invalid ESD measurement " + path.string() + ": " + error.what());
    }
}

std::string format_esd_table(const Measurement &m, const std::string &source) {
    const auto &[center, cosmology, sample, radial_bins] = m.metadata;
    std::ostringstream out;
    out << std::scientific << std::setprecision(16)
        << "# R_low R_high R_eff n_pair ESD_meas ESD_sigma ESD_cross ESD_cross_sigma\n"
        << "# Comoving R [Mpc/h], ESD [h M_sun/pc^2]\n"
        << "# center=" << bcg::center_name(center) << " groups=" << sample.n_lenses << " z_eff=" << sample.redshift
        << " mean_Mstar_Msun_over_h=" << sample.stellar_mass << " n_jack=" << sample.regions << '\n'
        << "# measurement=" << source << '\n'
        << "# FD restored: " << m.fd_restored << ", calibration: " << m.calibration.description() << '\n'
        << "# sources_after_cuts=" << m.sources_after_cuts << '\n';
    if (m.photometry_selected >= 0) out << "# photometry_select_ge_1=" << m.photometry_selected << '\n';
    if (m.calibration_sources >= 0)
        out << "# Self-calibration: Cartesian unique all-R union, sources=" << m.calibration_sources
            << ", FD interval=[" << -m.fd_limit << ',' << m.fd_limit << "], bins=" << m.fd_bins << '\n';
    for (const auto &bin : m.bins)
        out << bin.low << ' ' << bin.high << ' ' << bin.effective << ' ' << bin.count << ' '
            << bin.measured << ' ' << bin.sigma << ' ' << bin.cross << ' ' << bin.cross_sigma << '\n';
    return out.str();
}

void print_measurement_summary(const Measurement &m) {
    const auto &sample = m.metadata.sample;
    std::int64_t pairs = 0;
    for (const auto &bin : m.bins) {
        if (bin.count > std::numeric_limits<std::int64_t>::max() - pairs)
            throw std::runtime_error("Total lens–source pair count exceeds int64");
        pairs += bin.count;
    }
    std::cout << calibration_summary(m) << '\n'
              << "Sources after cuts: " << m.sources_after_cuts << '\n'
              << "Lens-source pairs: " << pairs << '\n'
              << "Combined " << bcg::center_name(m.metadata.center) << " stack: " << sample.n_lenses
              << " groups, mean z=" << sample.redshift;
    if (m.photometry_selected >= 0) std::cout << ", select>=1: " << m.photometry_selected;
    if (m.metadata.center == bcg::Center::bcg)
        std::cout << ", mean stellar mass [M_sun/h]=" << sample.stellar_mass;
    const auto range = std::minmax_element(m.groups_per_region.begin(), m.groups_per_region.end());
    std::cout << "\nJackknife regions: " << sample.regions << '/' << m.requested_regions
              << ", groups per region " << *range.first << " to " << *range.second << std::endl;
}
} // namespace bcg_esd
