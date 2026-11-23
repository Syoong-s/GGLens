#include "bcg_jackknife.hpp"
#include "bcg_pdf_sym.hpp"

namespace bcg {
// ==========================================
// Function: prepare_replica_pdf
// Method: Estimate equal-occupancy edges for all deletions from collective log-|G| histograms.
// ==========================================
ReplicaPdf prepare_replica_pdf(const PairStore &pairs, const std::vector<int> &labels,
                                int regions, int component, int bins, MPI_Comm comm) {
    ReplicaPdf pdf;
    pdf.edges.resize(regions + 1);
    pdf.moments.resize(regions + 1);
    pdf.kept.resize(regions + 1);
    double min_abs = std::numeric_limits<double>::infinity(), max_abs = 0;
    pairs.scan([&](const PairRecord &p) {
        const double value = component ? p.gx : p.gt, response = component ? p.dex : p.det;
        const double abs_value = std::abs(value);
        if (abs_value > 0) min_abs = std::min(min_abs, abs_value);
        max_abs = std::max(max_abs, abs_value);
        const int slot = 1 + labels[p.lens];
        ++pdf.kept[slot];
        pdf.moments[slot][0] += response * response;
        pdf.moments[slot][1] += response * value;
        pdf.moments[slot][2] += value * value;
    });
    MPI_Allreduce(MPI_IN_PLACE, &min_abs, 1, MPI_DOUBLE, MPI_MIN, comm);
    MPI_Allreduce(MPI_IN_PLACE, &max_abs, 1, MPI_DOUBLE, MPI_MAX, comm);
    MPI_Allreduce(MPI_IN_PLACE, pdf.kept.data(), regions + 1, MPI_INT64_T, MPI_SUM, comm);
    MPI_Allreduce(MPI_IN_PLACE, pdf.moments.data(), 3 * (regions + 1), MPI_DOUBLE, MPI_SUM, comm);
    for (int s = 1; s <= regions; ++s) {
        pdf.kept[0] += pdf.kept[s];
        for (int j = 0; j < 3; ++j) pdf.moments[0][j] += pdf.moments[s][j];
    }
    for (int s = 1; s <= regions; ++s) {
        pdf.kept[s] = pdf.kept[0] - pdf.kept[s];
        for (int j = 0; j < 3; ++j) pdf.moments[s][j] = pdf.moments[0][j] - pdf.moments[s][j];
    }

    constexpr int slices = 2048, stride = slices + 1; // Cell zero stores exact zeros.
    const double origin = max_abs > 0 ? std::log(min_abs) : 0;
    const double step = max_abs > min_abs ? (std::log(max_abs) - origin) / slices : 0;
    std::vector<std::int64_t> histogram(static_cast<std::size_t>(regions + 1) * stride);
    pairs.scan([&](const PairRecord &p) {
        const double value = std::abs(component ? p.gx : p.gt);
        int cell = 0;
        if (value > 0)
            cell = 1 + (step > 0 ? std::clamp(static_cast<int>((std::log(value) - origin) / step),
                                               0, slices - 1) : 0);
        ++histogram[(1 + labels[p.lens]) * stride + cell];
    });
    MPI_Allreduce(MPI_IN_PLACE, histogram.data(), static_cast<int>(histogram.size()),
                  MPI_INT64_T, MPI_SUM, comm);
    for (int s = 1; s <= regions; ++s)
        for (int cell = 0; cell < stride; ++cell) histogram[cell] += histogram[s * stride + cell];
    for (int s = 0; s <= regions; ++s) {
        if (pdf.kept[s] < bins + 2LL || !(pdf.moments[s][0] > 0)) continue;
        const auto *omitted = s > 0 ? histogram.data() + s * stride : nullptr;
        const auto count = [&](int cell) { return histogram[cell] - (omitted ? omitted[cell] : 0); };
        int cell = 0;
        std::int64_t before = 0;
        auto &edges = pdf.edges[s];
        edges.reserve(bins);
        for (int i = 1; i <= bins; ++i) {
            const double target = static_cast<double>(pdf.kept[s]) * i / (bins + 1);
            while (cell < slices && before + count(cell) < target) before += count(cell++);
            if (cell == 0) edges.push_back(0);
            else {
                const double fraction = std::clamp((target - before) / std::max<std::int64_t>(1, count(cell)),
                                                   0.0, 1.0);
                edges.push_back(std::min(max_abs, std::exp(origin + step * (cell - 1 + fraction))));
            }
        }
    }
    return pdf;
}
} // namespace bcg
