#include <math.h>
#include "Constants.hpp"
#include "ostream.hpp"

#include "mpi.h"
#include "fftPair/fftPair.hpp"

#include <fftw3-mpi.h>

/* we use the fourier transform of 
     f(i) = std::exp(-a*i*i)
   given by 
     F(w) = std::exp(-pi*pi*w*w/a)    
     
*/
int main(int argc, char **argv)
{
    size_t Npoints = 1000;
    int irank;
    int howmany = 1;

/*     ***********************get splitted arrays************************************************/
    std::ptrdiff_t alloc_local, local_n0, local_0_start;
    std::ptrdiff_t dimensions[] = {int(Npoints), 1};
    int required = MPI_THREAD_SINGLE;
    int provided;
    //ierr = MPI_Init ( &argc, &argv );
    auto ierr = MPI_Init_thread ( &argc, &argv, required, &provided );std::cout << "initialized.\n";
if ( ierr != 0 )
{
    std::cout << "\n";
    std::cout << "BONES - Fatal error!\n";
    std::cout << "MPI_Init returned ierr = " << ierr << "\n";
    return 0;
}    
fftw_mpi_init();    
std::cout << "fftw_mpi_init.\n";
    int rank = 2;
    int nproc;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &irank);
    std::cout << "Hello from rank " << irank << std::endl;

    //some care: this function does not work for 1d fft.
    alloc_local = fftw_mpi_local_size_many(rank, dimensions, howmany, 
                                   FFTW_MPI_DEFAULT_BLOCK, MPI_COMM_WORLD,
                                   &local_n0, &local_0_start);
std::cout << "after fftw_mpi_local_size_many.\n";

    std::cout << "local_n0: " << local_n0 << "local_0_start: " << local_0_start << std::endl;

    //auto data = fftw_alloc_complex(alloc_local);
/*     ***********************end of get splitted arrays************************************************/


    mdarray<std::complex<double>,2> Array_x({1,size_t(local_n0)});
    mdarray<std::complex<double>,2> Array_k({1,size_t(local_n0)});


    double a = 1./double(Npoints);
    for(int i=0; i<Npoints; i++){
        double x;
        if(i<=Npoints/2){
            x = i;
        }
        else{
            x = Npoints-i;
        }
        if( i>= local_0_start && i < local_n0+local_0_start ) {
            Array_x(0,i-int(local_0_start)) = std::exp(-a*std::pow(x, 2));
        }
    }

    std::cout << "Setting up Fourier Transform..\n";
{
    FourierTransform fftm(Array_x, Array_k, {int(Npoints),1});

    std::cout << "Doing Fourier transform...\n";
    fftm.fft(-1);
    std::cout << "DONE!\n";

    MPI_Barrier(MPI_COMM_WORLD);

    double DeltaW = 1./double(Npoints);
    std::ofstream os_rank;
    std::stringstream rankname;
    rankname << irank << ".txt";
    std::cout << rankname.str() << std::endl;
    os_rank.open(rankname.str());

    os_rank << "+------------+--------------------+---------------------+-------------------+\n";
    os_rank << "|    freq    | Numerical solution | Analytical solution |      Error(%)     |\n";
    os_rank << "+------------+--------------------+---------------------+-------------------+\n";

    auto&& NumericalSolution = fftm.get_Array_k();


    for(int i=0; i<Npoints; i++){
        if( i>= local_0_start && i < local_n0+local_0_start ) {
            double w;
            if(i<=Npoints/2){
                w = double(i)*DeltaW;
            }
            else{
                w = double(Npoints-i)*DeltaW;
            }
            auto AnalyticalSolution = std::sqrt(pi/a)*std::exp(-pi*pi*w*w/a)/std::sqrt(Npoints);
            os_rank << "|";
            os_rank << std::setw(7) << std::fixed << int(i);
            os_rank << "  |  ";
            os_rank << std::setw(6) << std::setprecision(2) <<  std::scientific << w;
            os_rank << "  ";
            os_rank << "|";
            os_rank << "  ";
            os_rank << std::setw(16) << std::setprecision(8) << std::scientific <<  NumericalSolution(0,i-int(local_0_start));
            os_rank << "  ";
            os_rank << "|";
            os_rank << "  ";
            os_rank << std::setw(16) << std::setprecision(8) << std::scientific << AnalyticalSolution;
            os_rank << "   ";
            os_rank << "|";
            os_rank << "  ";
            //os_rank << i-int(local_0_start) << std::endl;
            os_rank << std::setw(15) << std::setprecision(8) << std::scientific <<  100*abs(NumericalSolution(0,i-int(local_0_start))-AnalyticalSolution)/abs(AnalyticalSolution);
            os_rank << "  |" << std::endl;

            //if(AnalyticalSolution > 1.e-07 && 100*abs(NumericalSolution(0,i)-AnalyticalSolution)/abs(AnalyticalSolution)>5.){
            //   exit(1);
            //}
        }
    }
    os_rank.close();
    MPI_Barrier(MPI_COMM_WORLD);
}
    MPI_Finalize(); 

}
