#include <complex.h>
#include <math.h>

#include "Constants.hpp"
#include "ostream.hpp"
#include "mpi.h"
#include "fftPair/fftPair.hpp"
#include "MPIindex/MPIindex.hpp"
#include <fftw3-mpi.h>
#include <ios> 
/* we use the fourier transform of 
     f(i) = std::exp(-a*i*i)
   given by 
     F(w) = std::exp(-pi*pi*w*w/a)    
     
*/
const ptrdiff_t N0 = 200;
const ptrdiff_t N1 = 200;
const double a = 1./double(N0);
const double b = 1./double(N1);



int main(int argc, char **argv)
{
    int irank, nproc;
    int howmany = 1;

/*     ***********************get splitted arrays************************************************/
    std::ptrdiff_t alloc_local, local_n0, local_0_start;
    std::ptrdiff_t dimensions[] = {int(N0), int(N1), 1};
    std::array<size_t,3> dims = {int(N0), int(N1), 1};

#ifdef NEGF_MPI
    int required = MPI_THREAD_SINGLE;
    int provided;
    auto ierr = MPI_Init_thread ( &argc, &argv, required, &provided );
    std::cout << "initialized.\n";
    if ( ierr != 0 )
    {
        std::cout << "\n";
        std::cout << "BONES - Fatal error!\n";
        std::cout << "MPI_Init returned ierr = " << ierr << "\n";
        return 0;
    }    
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &irank);

    std::cout << " I am rank " << irank << ". total number = " << nproc << std::endl;
    fftw_mpi_init();    
#else 
    irank =0;
    nproc = 1;
#endif

    std::cout << "going to initialize mpindex.\n";
    MPIindex<3> mpindex(dims);

    std::cout << "mpindex initialized.\n";
/*     ***********************end of get splitted arrays************************************************/

    //allocate array with recommended size
    mdarray<std::complex<double>,2> Array_x( { 1, size_t( mpindex.get_RecommendedAllocate_fftw() ) } );//({1,size_t(N0)*size_t(N1)});
    mdarray<std::complex<double>,2> Array_k( { 1, size_t( mpindex.get_RecommendedAllocate_fftw() ) } );//({1,size_t(N0)*size_t(N1)});

    //filling Array_x
    //loop over local indices
    auto& LocalRange = mpindex.get_LocalRange();
    //#pragma omp parallel for
    for( int oneDindex_loc = 0; oneDindex_loc <= LocalRange.second - LocalRange.first; ++oneDindex_loc ) {
        auto aux_ = mpindex.loc1D_to_globnD(oneDindex_loc);
        //std::cout << oneDindex_loc << "/" << LocalRange.second << std::endl;
        double x, y;
        if( aux_[0] <= N0/2 )     x = aux_[0];
        else                      x = N0-aux_[0];
        if( aux_[1] <= N1/2 )     y = aux_[1];
        else                      y = N1-aux_[1];
  
        Array_x(0, oneDindex_loc) = std::exp(-a*std::pow(x, 2))*std::exp(-b*std::pow(y,2));
    }
    std::cout << "Rank " << irank <<  " Settingv up Fourier Transform..\n";
{
    FourierTransform fftm(Array_x, Array_k, {int(N0),int(N1),1});
    std::cout << "Rank " << irank <<  " Constructed Fourier Transform..\n";
    std::cout << "Doing Fourier transform...\n";
    fftm.fft(-1);
    std::cout << "DONE!\n";

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


    for( int oneDindex_loc = 0; oneDindex_loc <= LocalRange.second-LocalRange.first; ++oneDindex_loc ) {
        auto aux_ = mpindex.loc1D_to_globnD(oneDindex_loc);
        double wx, wy;
        if( aux_[0] <= N0/2 )     wx = double(aux_[0]) * DeltaWx;
        else                      wx = double(N0-aux_[0]) * DeltaWx;
        if(aux_[1] <= N1/2 )      wy = double(aux_[1]) * DeltaWy;
        else                      wy = double(N1-aux_[1]) * DeltaWy;
  
        auto AnalyticalSolution = std::sqrt(pi/b)*std::sqrt(pi/a)*std::exp(-pi*pi*wx*wx/a)/std::sqrt(N0)*std::exp(-pi*pi*wy*wy/b)/std::sqrt(N1);
        auto NumericalSolution  = Array_k( 0, oneDindex_loc );
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
        if(AnalyticalSolution > 1.e-07 && RelativeError > 5.){
            exit(1);
        }
    }
    os_rank.close();

}//end of scope of fftw
#ifdef NEGF_MPI
    MPI_Finalize(); 
#endif
}
