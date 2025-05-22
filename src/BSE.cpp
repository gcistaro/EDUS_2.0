
#include <ctime>
#include <iostream>
#include "InputVariables/simulation_parameters.hpp"
#include "BetheSalpeter/Hamiltonian.hpp"
#include "core/print_timing.hpp"
#include "initialize.hpp"
#include "core/projectdir.hpp"


#define variable(x)  (#x)
#include "core/githash.hpp"

#include <typeinfo>

int main(int argc, char *argv[])
{
    PROFILE_START("EDUS");
    if ( argc != 2 ) {
        throw std::runtime_error("To start the program: EDUS <input.json>");
        exit(0);
    }

    initialize();
    auto ctx = std::make_shared<Simulation_parameters>();
    ctx->import(std::string(argv[1]));
    BS_Hamiltonian bse; 
    bse.initialize(ctx, 0);
    finalize();

    PROFILE_STOP("EDUS");
#ifdef EDUS_MPI
    if( mpi::Communicator::world().rank() == 0 )
#endif
    {
        print_timing(1);
    }
}
