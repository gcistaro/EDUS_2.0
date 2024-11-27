//this function must be called at the beginning of any miniapp to allow for 
//initialization of EDUS important quantities, such as:
//- header printing with git infos
//- MPI_COMM_WORLD
//- fftw initialization

#include "omp.h"
#include "mkl.h"
#include <thread>
#include <iomanip>
#include "initialize.hpp"


#ifdef EDUS_MPI
mpi::Communicator kpool_comm;
mpi::Communicator band_comm;
int NumberKpools;
#endif
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
    mkl_set_num_threads(omp_get_max_threads());
#ifdef EDUS_MPI
    mpi::Communicator::initialize(MPI_THREAD_FUNNELED);
    NumberKpools = mpi::Communicator::world().size();
    fftw_mpi_init();  

    if( mpi::Communicator::world().rank() == 0 ) {
#endif 
#ifdef EDUS_FFTWTHREADS
    fftw_init_threads();
#endif
    print_header();
    std::cout << "**************************************************    PARALLELIZATION RECAP     **********************************************\n";
        std::cout << "*    OpenMP  threads:   *     ";
        std::cout << std::left << std::setw(95) << omp_get_max_threads() << "*\n";
        std::cout << "*    MKL  threads:      *     ";
        std::cout << std::left << std::setw(95) << mkl_get_max_threads() << "*\n";
#ifdef EDUS_FFTWTHREADS
    fftw_plan_with_nthreads(omp_get_max_threads());
        std::cout << "*    fftw3 threads:     *     ";
        std::cout << std::left << std::setw(95) << fftw_planner_nthreads() << "*\n";
        std::cout << "******************************************************************************************************************************\n";
#endif
#ifdef EDUS_MPI
    }
#endif 

#ifdef EDUS_MPI

    assert( mpi::Communicator::world().size()%NumberKpools == 0 );//for now i just implemented a rectangular MPI grid
    
    //create k point communicator
    kpool_comm.generate( mpi::Communicator::world(), mpi::Communicator::world().rank()/NumberKpools );

    //create band communicator
    band_comm.generate( mpi::Communicator::world(), mpi::Communicator::world().rank()%NumberKpools );

    //recap of mpi ranks -- to be done only for high verbosity
    if ( mpi::Communicator::world().rank() != 0 ) {
        auto aux = mpi::Communicator::world().rank();
        mpi::Communicator::world().send( &aux, 0 );
        aux = kpool_comm.rank();
        mpi::Communicator::world().send( &aux, 0 );
        aux = kpool_comm.size();
        mpi::Communicator::world().send( &aux, 0 );
        aux = band_comm.rank();
        mpi::Communicator::world().send( &aux, 0 );
        aux = band_comm.size();
        mpi::Communicator::world().send( &aux, 0 );
    }
    else {
        std::cout << "*    MPI world size:    *     ";
        std::cout << std::left << std::setw(95) << mpi::Communicator::world().size() << "*\n";
        //std::cout << "MPI parallelization recap. \n";
        //std::cout << "WORLD RANK: " << mpi::Communicator::world().rank() << "/" << mpi::Communicator::world().size();
        //std::cout << " KPOOL RANK: " << kpool_comm.rank() << "/" << kpool_comm.size();
        //std::cout << " BAND RANK: " << band_comm.rank()<< "/" << band_comm.size() << std::endl;

        for(int irank_ = 1; irank_ < mpi::Communicator::world().size(); ++irank_ ) {
            int world_comm_rank = 10, kpool_comm_rank, kpool_comm_size, band_comm_rank, band_comm_size;
            mpi::Communicator::world().receive( &world_comm_rank, irank_ );
            mpi::Communicator::world().receive( &kpool_comm_rank, irank_ );
            mpi::Communicator::world().receive( &kpool_comm_size, irank_ );
            mpi::Communicator::world().receive( &band_comm_rank, irank_ );
            mpi::Communicator::world().receive( &band_comm_size, irank_ );

            //std::cout << "WORLD RANK: " << world_comm_rank << "/" << mpi::Communicator::world().size();
            //std::cout << " KPOOL RANK: " << kpool_comm_rank << "/" << kpool_comm_size;
            //std::cout << " BAND RANK: " << band_comm_rank << "/" << band_comm_size << std::endl;
        }
    }
#endif
}



void finalize()
{  
#ifdef EDUS_MPI
    mpi::Communicator::finalize();
    if( mpi::Communicator::world().rank() == 0 ) {
        std::cout << "MPI Communicator finalized!\n";
        time_t now = time(0);
        char* dt = ctime(&now);
        std::cout << "Execution finished: " << dt << std::endl;
    }
#else
        time_t now = time(0);
        char* dt = ctime(&now);
        std::cout << "Execution finished: " << dt << std::endl;
#endif
}
