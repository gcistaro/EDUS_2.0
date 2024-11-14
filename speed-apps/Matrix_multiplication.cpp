#include "initialize.hpp"
#include "Geometry/Matrix.hpp"
#include "core/cmd_args/cmd_args.hpp"
#include "core/print_timing.hpp"


int main(int argn, char** argv)
{
    initialize();
    cmd_args args(argn, argv, 
                  {{"dim=", "{int} dimension of the matrix"},
                  {"niter=", "{int} number of matrix multiplication to average over"}}
                 );    


    PROFILE_START("Initialize")

    auto dim = args.value<int>("dim", 100);
    auto niter = args.value<int>("niter", 1000);

    Matrix<std::complex<double>> A(dim,dim); 
    Matrix<std::complex<double>> B(dim,dim); 
    Matrix<std::complex<double>> C(dim,dim); 

    A.fill(0.);
    B.fill(0.);
    C.fill(0.);

    PROFILE_STOP("Initialize")
    

    for(int icounter=0; icounter<niter; ++icounter) {
        PROFILE_START("Matrix Multiplication");
        Matrix_gemm(A, 1., B, C, 0.);  
        PROFILE_STOP("Matrix Multiplication");
    }

    print_timing(1);
}