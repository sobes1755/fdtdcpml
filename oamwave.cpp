#include "fdtdcpml.h"
#include "mdspanio.h"

#include <iostream>
#include <charconv>
#include <optional>
#include <cstring>
#include <fstream>
#include <numeric>

int
main(int argc, char *argv[])
{

    double x_min{0}, x_max{1};
    double y_min{0}, y_max{1};
    double z_min{0}, z_max{1};

    std::size_t p_int{257};
    std::size_t q_int{257};
    std::size_t r_int{257};

    if (argc != (1 + 9)) {
        std::cerr << "Usage: " << argv[0] << " xmin xmax pint ymin ymax qint zmin zmax rint\n";
        return EXIT_FAILURE;
    }

    std::from_chars(argv[1], argv[1] + std::strlen(argv[1]), x_min);
    std::from_chars(argv[2], argv[2] + std::strlen(argv[2]), x_max);
    std::from_chars(argv[3], argv[3] + std::strlen(argv[3]), p_int);
    std::from_chars(argv[4], argv[4] + std::strlen(argv[4]), y_min);
    std::from_chars(argv[5], argv[5] + std::strlen(argv[5]), y_max);
    std::from_chars(argv[6], argv[6] + std::strlen(argv[6]), q_int);
    std::from_chars(argv[7], argv[7] + std::strlen(argv[7]), z_min);
    std::from_chars(argv[8], argv[8] + std::strlen(argv[8]), z_max);
    std::from_chars(argv[9], argv[9] + std::strlen(argv[9]), r_int);

    std::optional<FDTDCPML> solver;

    try {
        solver.emplace(x_min, x_max, p_int, y_min, y_max, q_int, z_min, z_max, r_int);
    } catch (const std::exception& e) {
        std::cerr << "Memory allocation failed: " << e.what() << "\n";
        return EXIT_FAILURE;
    }

    // Solving and saving middle slices with Ex values...

    double x0 = std::midpoint(x_min, x_max);
    double y0 = std::midpoint(y_min, y_max);
    double z0 = std::midpoint(z_min, z_max);

    auto [p0, q0, r0] = solver->Ex_xyz2pqr({x0, y0, z0});    

    solver->init();

    for (size_t s = 0; s < 32; ++s) {

        solver->step();

        if (std::ofstream ofs{std::format("{}.{:d}.bin", "../out.bin/Ex.YZ", s), std::ios::binary}) {
            fdtdcpml::write_mdspan(ofs, solver->slice(FDTDCPML::ECOMP::Ex, FDTDCPML::EPQRS::P, p0));
        }
        if (std::ofstream ofs{std::format("{}.{:d}.bin", "../out.bin/Ex.XZ", s), std::ios::binary}) {
            fdtdcpml::write_mdspan(ofs, solver->slice(FDTDCPML::ECOMP::Ex, FDTDCPML::EPQRS::Q, q0));
        }
        if (std::ofstream ofs{std::format("{}.{:d}.bin", "../out.bin/Ex.XY", s), std::ios::binary}) {
            fdtdcpml::write_mdspan(ofs, solver->slice(FDTDCPML::ECOMP::Ex, FDTDCPML::EPQRS::R, r0));
        }

    }

    return EXIT_SUCCESS;

}
