//this function must be called at the beginning of any miniapp to allow for 
//initialization of NEGF important quantities, such as:
//- header printing with git infos
//- MPI_COMM_WORLD
//- fftw initialization
#include <sstream>
#include <fstream>

#include "core/print_header.hpp"
#include "core/mpi/Communicator.hpp"
#include "fftw3-mpi.h"


mpi::Communicator kpool_comm;
mpi::Communicator band_comm;
int NumberKpools = 2;



/*
    Parallelization using grid / example with 4 kpools, 12 ranks

    rows are kpools, columns are bandpools. inside the grid, I put the rank on the world communicator.

                  bandpool
                0         1         2
           --|--------------------------
    kpool  0 |  0         4         8        
           1 |  1         5         9 
           2 |  2         6        10
           3 |  3         7        11
*/
void initialize()
{  
#ifdef NEGF_MPI
    mpi::Communicator::initialize(MPI_THREAD_MULTIPLE);
    if( mpi::Communicator::world().rank() == 0 ) {
        std::cout << "MPI Communicator world size: " << mpi::Communicator::world().size() << "processors\n";
    }
    fftw_mpi_init();  

    assert( mpi::Communicator::world().size()%NumberKpools == 0 );//for now i just implemented a rectangular MPI grid
    
    //create k point communicator
    kpool_comm.generate( mpi::Communicator::world(), mpi::Communicator::world().rank()/NumberKpools );

    //create band communicator
    band_comm.generate( mpi::Communicator::world(), mpi::Communicator::world().rank()%NumberKpools );

    //recap of mpi ranks -- to be done only for high verbosity
    if ( mpi::Communicator::world().rank() != 0 ) {
        mpi::Communicator::world().send( mpi::Communicator::world().rank(), 0 );
        mpi::Communicator::world().send( kpool_comm.rank(), 0 );
        mpi::Communicator::world().send( kpool_comm.size(), 0 );
        mpi::Communicator::world().send( band_comm.rank(), 0 );
        mpi::Communicator::world().send( band_comm.size(), 0 );
    }
    else {
        print_header();
        std::cout << "MPI parallelization recap. \n";
        std::cout << "WORLD RANK: " << mpi::Communicator::world().rank() << "/" << mpi::Communicator::world().size();
        std::cout << " KPOOL RANK: " << kpool_comm.rank() << "/" << kpool_comm.size();
        std::cout << " BAND RANK: " << band_comm.rank()<< "/" << band_comm.size() << std::endl;

        for(int irank_ = 1; irank_ < mpi::Communicator::world().size(); ++irank_ ) {
            int world_comm_rank = 10, kpool_comm_rank, kpool_comm_size, band_comm_rank, band_comm_size;
            mpi::Communicator::world().receive( world_comm_rank, irank_ );
            mpi::Communicator::world().receive( kpool_comm_rank, irank_ );
            mpi::Communicator::world().receive( kpool_comm_size, irank_ );
            mpi::Communicator::world().receive( band_comm_rank, irank_ );
            mpi::Communicator::world().receive( band_comm_size, irank_ );

            std::cout << "WORLD RANK: " << world_comm_rank << "/" << mpi::Communicator::world().size();
            std::cout << " KPOOL RANK: " << kpool_comm_rank << "/" << kpool_comm_size;
            std::cout << " BAND RANK: " << band_comm_rank << "/" << band_comm_size << std::endl;
        }
    }
#else  
    print_header();
#endif
}

