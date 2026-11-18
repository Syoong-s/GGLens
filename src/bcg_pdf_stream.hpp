#pragma once
#include "bcg_pair_store.hpp"
#include "bcg_pdf_sym.hpp"

namespace bcg {
// ==========================================
// Function: pdf_sym_stream
// Method: Preserve PDF_SYM statistics while accepting a memory-bounded pair view.
// ==========================================
PdfResult pdf_sym_stream(const PairView &view, const std::vector<double> &initial, int nbin,
                        MPI_Comm comm, double center, double width, std::int64_t ntot_ref = -1);
} // namespace bcg
