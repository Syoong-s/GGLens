#include "bcg_dr9_rows.hpp"
#include <iostream>

namespace bcg {
namespace {
int first_row_width(const std::string &path, const CatalogTextFormat &format) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("Cannot open catalog: " + path);
    if (format.header) in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::string line;
    while (std::getline(in, line)) {
        if (!is_catalog_data_row(line)) continue;
        std::string_view remaining = line;
        std::size_t width = 0;
        while (take_catalog_field(remaining, format.delimiter).data()) ++width;
        if (width > static_cast<std::size_t>(std::numeric_limits<int>::max()))
            throw std::runtime_error("Too many catalog columns: " + path);
        return static_cast<int>(width);
    }
    if (!in.eof() || in.bad()) throw std::runtime_error("Read failed for catalog: " + path);
    throw std::runtime_error("No data row in first catalog: " + path);
}

int resolve_column(int column, int width, const char *key, bool allow_zero = false) {
    if (column == 0 && allow_zero) return 0;
    const int resolved = column < 0 ? width + column + 1 : column;
    if (resolved < 1 || resolved > width) {
        throw std::runtime_error(std::string(key) + " = " + std::to_string(column) +
                                 " is outside the catalog's " + std::to_string(width) + " columns");
    }
    return resolved;
}
} // namespace

LensColumns read_lens_columns(const std::vector<CatalogTriple> &triples,
                              LensColumns columns, Center center, MPI_Comm comm) {
    int rank = 0;
    MPI_Comm_rank(comm, &rank);
    std::array<int, 3> widths{};
    auto &[group, galaxy, member] = columns;
    if (rank == 0) {
        const auto &first = triples.front();
        widths = {first_row_width(first.group, group.text), first_row_width(first.galaxy, galaxy.text),
                  first_row_width(first.member, member.text)};
        std::cout << "First catalog columns: Group " << widths[0] << ", Galaxy " << widths[1]
                  << ", Member " << widths[2] << std::endl;
    }
    MPI_Bcast(widths.data(), static_cast<int>(widths.size()), MPI_INT, 0, comm);
    const auto [group_width, galaxy_width, member_width] = widths;
    group.id = resolve_column(group.id, group_width, "group_columns.id", true);
    group.redshift = resolve_column(group.redshift, group_width, "group_columns.redshift");
    group.halo_mass = resolve_column(group.halo_mass, group_width, "group_columns.log-halo-mass");
    galaxy.id = resolve_column(galaxy.id, galaxy_width, "galaxy_columns.id", true);
    galaxy.ra = resolve_column(galaxy.ra, galaxy_width, "galaxy_columns.ra");
    galaxy.dec = resolve_column(galaxy.dec, galaxy_width, "galaxy_columns.dec");
    member.galaxy_id = resolve_column(member.galaxy_id, member_width, "member_columns.galaxy-id");
    member.group_id = resolve_column(member.group_id, member_width, "member_columns.group-id");
    galaxy.luminosity = center == Center::lwc ?
        resolve_column(galaxy.luminosity, galaxy_width, "galaxy_columns.luminosity") : 0;
    if (center == Center::bcg) {
        galaxy.stellar_mass = resolve_column(galaxy.stellar_mass, galaxy_width, "galaxy_columns.log-stellar-mass");
        galaxy.photometry_select = resolve_column(galaxy.photometry_select, galaxy_width,
                                                  "galaxy_columns.photometry-select", true);
        member.brightest = resolve_column(member.brightest, member_width, "member_columns.brightest-flag");
    } else {
        galaxy.stellar_mass = galaxy.photometry_select = member.brightest = 0;
    }
    return columns;
}

std::int64_t preceding_catalog_rows(std::int64_t local_rows, MPI_Comm comm) {
    int rank = 0;
    MPI_Comm_rank(comm, &rank);
    std::int64_t offset = 0;
    MPI_Exscan(&local_rows, &offset, 1, MPI_INT64_T, MPI_SUM, comm);
    return rank == 0 ? 0 : offset;
}

std::int64_t galaxy_row_offset(const std::string &path, const GalaxyColumns &columns, MPI_Comm comm) {
    if (columns.id != 0) return 0;
    int rank = 0, ranks = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &ranks);
    if (ranks == 1) return 0;

    // Galaxy matches need the global ID before parsing any fields.
    std::int64_t rows = 0;
    if (rank < ranks - 1) {
        Dr9Rows input(path, columns.text.header, comm);
        std::string line;
        while (input.next(line)) {}
        rows = input.row_number();
    }
    return preceding_catalog_rows(rows, comm);
}
} // namespace bcg
