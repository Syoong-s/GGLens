#pragma once

#include "bcg_calibration.hpp"
#include "bcg_physics.hpp"
#include "bcg_stream.hpp"

namespace bcg {
struct SourceMoments {
    double g1, g2, de, h1, h2;
};
static_assert(sizeof(SourceMoments) == 40, "A matched source stores five double moments");

struct LensGeometry {
    double chi = 0, chord = 0, z = 0;
    std::array<double, 3> xyz{};
};

struct PreparedSource {
    SourceMoments moments;
    double zp, chi;
    std::array<double, 3> xyz;
};

// Share one rank-local table across all R bins. A pair indexes only within its
// source chunk, so a rank may process more than 2^32 sources in total.
class SourceMomentStore {
#pragma pack(push, 4)
    struct DiskSource {
        double g1, g2, de, h1, h2;
        double ra, dec, zp;
        std::uint32_t order;
    };
    struct RawSource {
        DiskSource source;
        double gf1, gf2;
    };
#pragma pack(pop)
    static_assert(sizeof(DiskSource) == 68, "Store eight doubles and a source index");
    static_assert(sizeof(RawSource) == 84, "Retain two additional doubles until self-calibration");
    struct Block {
        std::uint64_t start;
        std::uint32_t count, slots;
    };

    std::string file_;
    std::size_t io_bytes_;
    std::unique_ptr<BinaryWriter<DiskSource>> writer_;
    std::unique_ptr<BinaryWriter<RawSource>> raw_writer_;
    std::uint64_t rows_ = 0, block_start_ = 0;
    std::vector<Block> blocks_;
    std::ifstream reader_;
    std::size_t loaded_ = std::numeric_limits<std::size_t>::max();
    std::vector<PreparedSource> memory_;
    Cosmology cosmology_;
    std::vector<LensGeometry> lenses_;

public:
    SourceMomentStore(const std::string &file, std::size_t io_bytes, bool self_calibrate,
                      const Cosmology &cosmology);
    ~SourceMomentStore();
    std::uint32_t block_id() const;
    void append(const Source &source, std::uint32_t order);
    void finish_block(std::size_t slots);
    void finish(std::vector<LensGeometry> lenses);
    void calibrate(const ShearCalibration &calibration);
    const std::vector<PreparedSource> &load_block(std::uint32_t id);
    const LensGeometry &lens(std::uint32_t i) const { return lenses_.at(i); }
};
} // namespace bcg
