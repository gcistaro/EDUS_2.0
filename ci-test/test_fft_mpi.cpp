#include <complex.h>
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
const ptrdiff_t N0 = 300;
const ptrdiff_t N1 = 300;
const double a = 1./double(N0);
const double b = 1./double(N1);
//std::shared_ptr<fftw_complex> my_function(const int& i, const int& j)
//{
//        double x;
//        if(i<=N0/2){
//            x = i;
//        }
//        else{
//            x = N0-i;
//        }
//        fftw_complex f;
//        f = std::exp(-a*i*i);
//        return std::make_shared<fftw_complex>(f);
//}

//int main(int argc, char **argv)
//{
//    fftw_plan plan;
//    fftw_complex *data;
//    ptrdiff_t alloc_local, local_n0, local_0_start, i, j;
//
//    MPI_Init(&argc, &argv);
//    fftw_mpi_init();
//
//    /* get local data size and allocate */
//    alloc_local = fftw_mpi_local_size_2d(N0, N1, MPI_COMM_WORLD,
//                                         &local_n0, &local_0_start);
//std::cout << "local_n0: " << local_n0 << "local_0_start: " << local_0_start << " alloc_local: " << alloc_local << std::endl;
//    data = fftw_alloc_complex(alloc_local);
//
//    /* create plan for in-place forward DFT */
//    plan = fftw_mpi_plan_dft_2d(N0, N1, data, data, MPI_COMM_WORLD,
//                                FFTW_FORWARD, FFTW_ESTIMATE);    
//
//    /* initialize data to some function my_function(x,y) */
//    //for (i = 0; i < local_n0; ++i) for (j = 0; j < N1; ++j)
//    //   data[i*N1 + j] = (*(my_function(local_0_start + i, j)));
//
//    /* compute transforms, in-place, as many times as desired */
//    fftw_execute(plan);
//
//    fftw_destroy_plan(plan);
//
//    MPI_Finalize();
//}

int main(int argc, char **argv)
{
    int irank;
    int howmany = 1;

/*     ***********************get splitted arrays************************************************/
    std::ptrdiff_t alloc_local, local_n0, local_0_start;
    std::ptrdiff_t dimensions[] = {int(N0), int(N1)};
    int required = MPI_THREAD_SINGLE;
    int provided;
    //ierr = MPI_Init ( &argc, &argv );
    auto ierr = MPI_Init_thread ( &argc, &argv, required, &provided );
    std::cout << "initialized.\n";
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

    std::cout << "local_n0: " << local_n0 << "local_0_start: " << local_0_start << " alloc_local: " << alloc_local << std::endl;

    //auto data = fftw_alloc_complex(alloc_local);
/*     ***********************end of get splitted arrays************************************************/


    mdarray<std::complex<double>,2> Array_x({1,size_t(alloc_local)});//({1,size_t(N0)*size_t(N1)});
    mdarray<std::complex<double>,2> Array_k({1,size_t(alloc_local)});//({1,size_t(N0)*size_t(N1)});

  //filling Array_x
    for(int i=0; i<N0; i++){
        double x;
        if(i<=N0/2) {
          x = i;
        }
        else {
          x = N0-i;
        }
        if( i>= local_0_start && i < local_n0+local_0_start ) {          
            for(int j=0; j<N1; j++) {
                double y;
                if(j<=N1/2) {
                    y = j;
                }
                else {
                    y = N1-j;
                }
                
                Array_x(0,j+N1*int(i-int(local_0_start))) = std::exp(-a*std::pow(x, 2))*std::exp(-b*std::pow(y,2));
            }
        }
    }

    std::cout << "Rank " << irank <<  " Setting up Fourier Transform..\n";
{
    FourierTransform fftm(Array_x, Array_k, {int(N0),int(N1)});
    std::cout << "Rank " << irank <<  " Constructed Fourier Transform..\n";
    MPI_Barrier(MPI_COMM_WORLD);

    std::cout << "Doing Fourier transform...\n";
    fftm.fft(-1);
    std::cout << "DONE!\n";

    MPI_Barrier(MPI_COMM_WORLD);

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

    auto&& NumericalSolution = fftm.get_Array_k();


    for(int i=0; i<N0; i++){
        if( i>= local_0_start && i < local_n0+local_0_start ) {
            double wx;
            if(i<=N0/2){
                wx = double(i)*DeltaWx;
            }
            else{
                wx = double(N0-i)*DeltaWx;
            }

            for( int j=0; j < N1; j++ ) {
                double wy;
                if(j<=N1/2){
                    wy = double(j)*DeltaWy;
                }
                else{
                    wy = double(N1-j)*DeltaWy;
                }
                auto AnalyticalSolution = std::sqrt(pi/b)*std::sqrt(pi/a)*std::exp(-pi*pi*wx*wx/a)/std::sqrt(N0)*std::exp(-pi*pi*wy*wy/b)/std::sqrt(N1);
                os_rank << "|";
                os_rank << std::setw(7) << std::fixed << int(i);
                os_rank << "  |  ";
                os_rank << std::setw(6) << std::setprecision(2) <<  std::scientific << wx;
                os_rank << "  ";
                os_rank << "|";
                os_rank << "  ";
                os_rank << std::setw(16) << std::setprecision(8) << std::scientific <<  NumericalSolution(0,j+N1*int(i-int(local_0_start))) ;
                os_rank << "  ";
                os_rank << "|";
                os_rank << "  ";
                os_rank << std::setw(16) << std::setprecision(8) << std::scientific << AnalyticalSolution;
                os_rank << "   ";
                os_rank << "|";
                os_rank << "  ";
                //os_rank << i-int(local_0_start) << std::endl;
                os_rank << std::setw(15) << std::setprecision(8) << std::scientific <<  100*abs(NumericalSolution(0,j+N1*int(i-int(local_0_start)))-AnalyticalSolution )/abs(AnalyticalSolution);
                os_rank << "  |" << std::endl;
    
                if(AnalyticalSolution > 1.e-07 && 100*abs(NumericalSolution(0,j+N1*int(i-int(local_0_start)))-AnalyticalSolution )/abs(AnalyticalSolution)>5.){
                   std::cout << AnalyticalSolution << " " << 100*abs(NumericalSolution(0,i)-AnalyticalSolution)/abs(AnalyticalSolution) << std::endl;
                   exit(1);
                }
            }
        }
    }
    os_rank.close();
    MPI_Barrier(MPI_COMM_WORLD);

}//end of scope of fftw

    MPI_Finalize(); 

}
