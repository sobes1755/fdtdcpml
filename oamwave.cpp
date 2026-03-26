#include "fdtdcpml.h"
#include "mdio.h"

#include <iostream>
#include <charconv>
#include <optional>
#include <cstring>
#include <fstream>

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

    //

    std::cout << "x_step = " << solver->x_step() << std::endl;
    std::cout << "y_step = " << solver->y_step() << std::endl;
    std::cout << "z_step = " << solver->z_step() << std::endl;
    std::cout << "t_step = " << solver->t_step() << std::endl;

    // FDTDCPML...

    solver->init();

    for (size_t s = 0; s < 1024; ++s) {

//        std::cout << "s = " << s << ", t = " << solver->t() << std::endl;

        solver->step();

//        if (s % 16 == 0) {

            std::ofstream ofs(std::format("{}.{:d}.bin", "../out.bin/Ex", s), std::ios::binary);

            if (ofs.is_open()) {

//                auto sliceRaw = std::submdspan(solver->Ez,
//                    std::strided_slice{.offset = 0, .extent = solver->sizeX, .stride = 1},
//                    std::strided_slice{.offset = 0, .extent = solver->sizeY, .stride = 1},
//                    rickerEzZ);

                auto sliceRaw = std::submdspan(solver->Ex(),
                    std::full_extent,
                    std::full_extent,
                    61);

                auto sliceScaled = std::mdspan(sliceRaw.data_handle(),
                    sliceRaw.mapping(), scaling_accessor<double>{solver->x_step()});

                writeSlice(ofs, sliceScaled);

                ofs.close();

//            }

        }

    }

    return EXIT_SUCCESS;

}
