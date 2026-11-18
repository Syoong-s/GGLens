#pragma once
#include "selfcalibration.hpp"

namespace bcg_esd {
void render_plots(const std::vector<Result> &results, const Config &cfg,
                  const SelfCalibrationDiagnostics &diagnostics);
} // namespace bcg_esd
