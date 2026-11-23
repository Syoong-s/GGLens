#include "bcg_esd.hpp"
#include <iostream>

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    try {
        if (argc != 2) throw std::runtime_error("Expected one INI path");
        const auto cfg = bcg_esd::load_config(argv[1], MPI_COMM_WORLD);
        const auto lenses = bcg_esd::load_lenses(cfg, MPI_COMM_WORLD);
        if (rank == 0) {
            std::cout << std::setprecision(17);
            for (const auto &lens : lenses)
                std::cout << "lens " << lens.catalog_id << ' ' << lens.id << ' ' << lens.ra << ' '
                          << lens.dec << ' ' << lens.z << ' ' << lens.log_mh << ' ' << lens.galaxy_id
                          << ' ' << lens.stellar_mass << ' ' << lens.photometry_select << '\n';
        }
        MPI_Finalize();
    } catch (const std::exception &e) {
        std::cerr << "Rank " << rank << " error: " << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
}
