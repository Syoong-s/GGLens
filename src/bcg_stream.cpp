#include "bcg_stream.hpp"

namespace bcg {
// ==========================================
// Function: ScratchSpace constructor
// Method: Broadcast a unique run name; create rank directories without overwrites.
// ==========================================
ScratchSpace::ScratchSpace(const std::string &base, MPI_Comm comm) {
    int rank = 0; MPI_Comm_rank(comm, &rank);
    std::string token;
    if (!rank) token = "run_" + std::to_string(
        std::chrono::system_clock::now().time_since_epoch().count()) + "_" + std::to_string(getpid());
    int n = static_cast<int>(token.size()); MPI_Bcast(&n, 1, MPI_INT, 0, comm);
    token.resize(n); MPI_Bcast(token.data(), n, MPI_CHAR, 0, comm);
    root_ = fs::absolute(base) / token;
    fs::create_directories(root_);
    rank_dir_ = root_ / ("rank_" + std::to_string(rank));
    if (!fs::create_directory(rank_dir_))
        throw std::runtime_error("Scratch rank directory already exists: " + rank_dir_.string());
    owned_ = true;
}

// ==========================================
// Function: ScratchSpace destructor
// Method: Best-effort cleanup of only this rank's owned directory during unwinding.
// ==========================================
ScratchSpace::~ScratchSpace() {
    if (owned_) { std::error_code error; fs::remove_all(rank_dir_, error); }
}

// ==========================================
// Function: ScratchSpace::path
// Method: Return a safe, unique scratch filename inside the owned directory.
// ==========================================
std::string ScratchSpace::path(const std::string &name) const {
    if (name.empty() || name.find('/') != std::string::npos || name == "." || name == "..")
        throw std::runtime_error("Invalid internal scratch filename");
    return (rank_dir_ / name).string();
}

// ==========================================
// Function: ScratchSpace::cleanup
// Method: Remove rank data, synchronize completion, and remove the empty run root.
// ==========================================
void ScratchSpace::cleanup(MPI_Comm comm) {
    fs::remove_all(rank_dir_); owned_ = false;
    MPI_Barrier(comm);
    int rank = 0; MPI_Comm_rank(comm, &rank);
    if (!rank) fs::remove(root_);
}
} // namespace bcg
