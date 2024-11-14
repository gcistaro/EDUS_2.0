#include "initialize.hpp"
#include "Operator/BlockMatrix.hpp"
#include "core/cmd_args/cmd_args.hpp"
#include "core/print_timing.hpp"


int main(int argn, char** argv)
{
    initialize();
    cmd_args args(argn, argv, 
                  {{"dim=", "{int} dimension of the matrix"},
                  {"niter=", "{int} number of matrix multiplication to average over"},
                  {"nblocks=", "(int) number of blocks of matrices to multiply"}}
                 );    


    PROFILE_START("Initialize")

    auto dim = args.value<int>("dim", 2);
    auto niter = args.value<int>("niter", 10000);
    auto nblocks = args.value<int>("nblocks", 1600);

    BlockMatrix<std::complex<double>> A(Space::k, nblocks,dim,dim); 
    BlockMatrix<std::complex<double>> B(Space::k, nblocks,dim,dim); 
    BlockMatrix<std::complex<double>> C(Space::k, nblocks,dim,dim); 

    A.fill(0.);
    B.fill(0.);
    C.fill(0.);

    PROFILE_STOP("Initialize")
    

    for(int icounter=0; icounter<niter; ++icounter) {
        PROFILE_START("BlockMatrix Multiplication");
        multiply(A, 1.+im*0., B, C, 0.*im);  
        PROFILE_STOP("BlockMatrix Multiplication");
    }

    print_timing(1);
}