#pragma once
#include "bcg_esd.hpp"

namespace bcg_esd {
struct FitChain {
    std::vector<double> samples, diagnostics;
};

// ==========================================
// Class: FitBackend
// Method: Keep embedded Python and the in-memory fit problem scoped to each sampling rank.
// ==========================================
class FitBackend {
    struct Runtime;
    std::unique_ptr<Runtime> runtime_;
public:
    explicit FitBackend(const Config &cfg);
    ~FitBackend();
    std::vector<double> correlation(const FitInput::Metadata &metadata);
    void setup(const FitInput &input, const std::vector<double> &corr);
    FitChain sample(const Config &cfg, int ensemble);
    void finish(const FitChain &chain, const Config &cfg, const FitInput &input);
};
} // namespace bcg_esd
