//this function must be called at the beginning of any miniapp to allow for 
//initialization of NEGF important quantities, such as:
//- header printing with git infos
//- MPI_COMM_WORLD
//- fftw initialization
#include "core/print_header.hpp"
#include "core/mpi/Communicator.hpp"
#include "fftw3-mpi.h"

void initialize()
{
    print_header();
#ifdef NEGF_MPI
    mpi::Communicator::initialize(MPI_THREAD_MULTIPLE);
    std::cout << " I am rank " << mpi::Communicator::world().rank();
    std::cout << ". total number = " << mpi::Communicator::world().size() << std::endl;
    fftw_mpi_init();  
#endif  
}

