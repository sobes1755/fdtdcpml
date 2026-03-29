#include "fdtdcpml.h"
#include "mdspanio.h"

#include <iostream>
#include <charconv>
#include <optional>
#include <cstring>
#include <fstream>
#include <numeric>
#include <print>

int
main(int argc, char *argv[])
{

    double x_min{0}, x_max{1};
    double y_min{0}, y_max{1};
    double z_min{0}, z_max{1};
    double t_min{0}, t_max{1};

    std::size_t p_int{257};
    std::size_t q_int{257};
    std::size_t r_int{257};
    std::size_t s_int{257};

    if (argc != (1 + 12)) {
        std::cerr << "Usage: " << argv[0] << " xmin xmax pint ymin ymax qint zmin zmax rint tmin tmax sint\n";
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
    std::from_chars(argv[7], argv[10] + std::strlen(argv[10]), t_min);
    std::from_chars(argv[8], argv[11] + std::strlen(argv[11]), t_max);
    std::from_chars(argv[9], argv[12] + std::strlen(argv[12]), s_int);

    std::optional<FDTDCPML> solver;

    try {
        solver.emplace(x_min, x_max, p_int, y_min, y_max, q_int, z_min, z_max, r_int, t_min, t_max, s_int);
    } catch (const std::exception& e) {
        std::cerr << "Memory allocation failed: " << e.what() << "\n";
        return EXIT_FAILURE;
    }

    // Solving and saving middle slices with Ex values...

    double x0 = std::midpoint(x_min, x_max);
    double y0 = std::midpoint(y_min, y_max);
    double z0 = std::midpoint(z_min, z_max);

    auto [p0, q0, r0] = solver->Ex_xyz2pqr({x0, y0, z0});

    std::print("x0 = {}, y0 = {}, z0 = {}.\n", x0, y0, z0);
    std::print("p0 = {}, q0 = {}, r0 = {}.\n", p0, q0, r0);

    solver->init();

    int frame = 0;

    auto print_cout = [&](){ std::print("t = {}, frame = {}.\n", solver->t(), frame); };

    auto print_ofs_binary = [&](){  
       if (std::ofstream ofs{std::format("{}.{:d}.bin", "../out.bin/Ex.YZ", frame), std::ios::binary}) {
           fdtdcpml::write_mdspan_binary(ofs, solver->slice(FDTDCPML::ECOMP::Ex, FDTDCPML::EPQRS::P, p0));
       }
       if (std::ofstream ofs{std::format("{}.{:d}.bin", "../out.bin/Ex.XZ", frame), std::ios::binary}) {
           fdtdcpml::write_mdspan_binary(ofs, solver->slice(FDTDCPML::ECOMP::Ex, FDTDCPML::EPQRS::Q, q0));
       }
       if (std::ofstream ofs{std::format("{}.{:d}.bin", "../out.bin/Ex.XY", frame), std::ios::binary}) {
           fdtdcpml::write_mdspan_binary(ofs, solver->slice(FDTDCPML::ECOMP::Ex, FDTDCPML::EPQRS::R, r0));
       }
    };

/*

    auto print_ofs_text = [&]() {
        if (std::ofstream ofs{std::format("{}.{:d}.txt", "../out.bin/Ex.YZ", frame)}) {
            fdtdcpml::write_mdspan_text(ofs, solver->slice(FDTDCPML::ECOMP::Ex, FDTDCPML::EPQRS::P, p0), std::format("x = {}, p = {}", x0, p0));
        }
        if (std::ofstream ofs{std::format("{}.{:d}.txt", "../out.bin/Ex.XZ", frame)}) {
            fdtdcpml::write_mdspan_text(ofs, solver->slice(FDTDCPML::ECOMP::Ex, FDTDCPML::EPQRS::Q, q0), std::format("y = {}, q = {}", y0, q0));
        }
        if (std::ofstream ofs{std::format("{}.{:d}.txt", "../out.bin/Ex.XY", frame)}) {
            fdtdcpml::write_mdspan_text(ofs, solver->slice(FDTDCPML::ECOMP::Ex, FDTDCPML::EPQRS::R, r0), std::format("z = {}, r = {}", z0, r0));
        }
    };

*/

    print_cout();
    print_ofs_binary();

    while (solver->t() <= t_max) {

        solver->step();

        ++frame;

        print_cout();
        print_ofs_binary();

    }

    return EXIT_SUCCESS;

}
