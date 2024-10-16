#include <complex.h>
#include <math.h>

#include "Constants.hpp"
#include "ostream.hpp"
#include "mpi.h"
#include "fftPair/fftPair.hpp"
#include "MPIindex/MPIindex.hpp"
#include "initialize.hpp"
#include <fftw3-mpi.h>
#include <ios> 
/* we use the fourier transform of 
     f(i) = std::exp(-a*i*i)
   given by 
     F(w) = std::exp(-pi*pi*w*w/a)    
     
*/
const ptrdiff_t N0 = 100;
const ptrdiff_t N1 = 100;
const double a = 1./double(N0);
const double b = 1./double(N1);



int main(int argc, char **argv)
{
    initialize();
    int irank, nproc;
    int howmany = 1;
/*     ***********************get splitted arrays************************************************/
    std::array<int,3> dims = {N0, N1, 1};

#ifdef NEGF_MPI
    irank = mpi::Communicator::world().rank();
    nproc = mpi::Communicator::world().size();
#else 
    irank =0;
    nproc = 1;
#endif

    std::cout << "going to initialize mpindex.\n";
    MPIindex<3> mpindex(dims);

    std::cout << "mpindex initialized.\n";
/*     ***********************end of get splitted arrays************************************************/

    //allocate array with recommended size
#ifdef NEGF_MPI
    mdarray<std::complex<double>,2> Array_x( { mpindex.get_RecommendedAllocate_fftw(), howmany  } );//({1,int(N0)*int(N1)});
    mdarray<std::complex<double>,2> Array_k( { mpindex.get_RecommendedAllocate_fftw(), howmany  } );//({1,int(N0)*int(N1)});
#else
    mdarray<std::complex<double>,2> Array_x( { howmany, int(N0)*int(N1)});
    mdarray<std::complex<double>,2> Array_k( { howmany, int(N0)*int(N1)});
#endif 

    //filling Array_x
    //loop over local indices
    #pragma omp parallel for
    for( int oneDindex_loc = 0; oneDindex_loc < mpindex.nlocal; ++oneDindex_loc ) {
        auto aux_ = mpindex.loc1D_to_globnD(oneDindex_loc);
        double x, y;
        if( aux_[0] <= N0/2 )     x = aux_[0];
        else                      x = N0-aux_[0];
        if( aux_[1] <= N1/2 )     y = aux_[1];
        else                      y = N1-aux_[1];
#ifdef NEGF_MPI
        Array_x(oneDindex_loc, 0) =
#else
        Array_x(0, oneDindex_loc) =
#endif 
        std::exp(-a*std::pow(x, 2))*std::exp(-b*std::pow(y,2));
    }
    std::cout << "Rank " << irank <<  " Settingv up Fourier Transform..\n";
{
    FourierTransform fftm(Array_x, Array_k, {int(N0),int(N1),1});
    std::cout << "Rank " << irank <<  " Constructed Fourier Transform..\n";
    std::cout << "Rank " << irank <<  "Doing Fourier transform...\n";
    fftm.fft(-1);
    std::cout << "Rank " << irank <<  "DONE Fourier transform!\n";

#ifdef NEGF_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    double DeltaWx = 1./double(N0);
    double DeltaWy = 1./double(N1);
    std::ofstream os_rank;
    std::stringstream rankname;
    rankname << irank << ".txt";
    std::cout << rankname.str() << std::endl;
    os_rank.open(rankname.str());

    os_rank << "+------------+--------------------+---------------------+-------------------+\n";
    os_rank << "|    freq    | Numerical solution | Analytical solution |      Error(%)     |\n";
    os_rank << "+------------+--------------------+---------------------+-------------------+\n";


    for( int oneDindex_loc = 0; oneDindex_loc < mpindex.nlocal; ++oneDindex_loc ) {
        auto aux_ = mpindex.loc1D_to_globnD(oneDindex_loc);
        double wx, wy;
        if( aux_[0] <= N0/2 )     wx = double(aux_[0]) * DeltaWx;
        else                      wx = double(N0-aux_[0]) * DeltaWx;
        if(aux_[1] <= N1/2 )      wy = double(aux_[1]) * DeltaWy;
        else                      wy = double(N1-aux_[1]) * DeltaWy;
  
        auto AnalyticalSolution = std::sqrt(pi/b)*std::sqrt(pi/a)*std::exp(-pi*pi*wx*wx/a)/std::sqrt(N0)*std::exp(-pi*pi*wy*wy/b)/std::sqrt(N1);
#ifdef NEGF_MPI
        auto NumericalSolution  = Array_k( oneDindex_loc, 0 );
#else
        auto NumericalSolution  = Array_k( 0, oneDindex_loc );
#endif
        auto RelativeError = 100*abs(NumericalSolution - AnalyticalSolution )/abs(AnalyticalSolution);
        //print all infos
        os_rank << "  |  ";
        os_rank << std::setw(6) << std::setprecision(2) <<  std::scientific << wx;
        os_rank << "  ";
        os_rank << "|";
        os_rank << "  ";
        os_rank << std::setw(6) << std::setprecision(2) <<  std::scientific << wy;
        os_rank << "  ";
        os_rank << "|";
        os_rank << "  ";
        os_rank << std::setw(16) << std::setprecision(8) << std::scientific <<  NumericalSolution;
        os_rank << "  ";
        os_rank << "|";
        os_rank << "  ";
        os_rank << std::setw(16) << std::setprecision(8) << std::scientific << AnalyticalSolution;
        os_rank << "   ";
        os_rank << "|";
        os_rank << "  ";
        os_rank << std::setw(15) << std::setprecision(8) << std::scientific << RelativeError;
        os_rank << std::endl;
        if(AnalyticalSolution > 1.e-07 && RelativeError > 0.5){
            exit(1);
        }
    }
    os_rank.close();

}//end of scope of fftw
#ifdef NEGF_MPI
    MPI_Finalize(); 
#endif
}
