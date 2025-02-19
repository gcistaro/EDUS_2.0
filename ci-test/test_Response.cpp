#include <ctime>
#include <iostream>
#include "InputVariables/simulation_parameters.hpp"
#include "Screened_potential/Screened_potential.hpp"
#include "core/print_timing.hpp"
#include "initialize.hpp"
#include "core/projectdir.hpp"

int main(int argc, char *argv[])
{
    PROFILE_START("EDUS");
    if ( argc != 2 ) {
        throw std::runtime_error("To start the program: EDUS <input.json>");
        exit(0);
    }
std::cout << "1"<<std::endl;
    initialize();

    auto ctx = std::make_shared<Simulation_parameters>();
    ctx->import(std::string(argv[1]));
    ScreenedPotential screened(ctx);

    std::cout << "omp_get_num_threads: " << omp_get_num_threads() << std::endl;
    screened.initialize();
    screened.Calculate();

    finalize();

    PROFILE_STOP("EDUS");

    print_timing(1);
}