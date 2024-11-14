#include "initialize.hpp"
#include "Operator/BlockMatrix.hpp"
#include "core/cmd_args/cmd_args.hpp"
#include "fftPair/fftPair.hpp"
#include "core/print_timing.hpp"


int main(int argn, char** argv)
{
    initialize();
    cmd_args args(argn, argv, 
                  {{"dim=", "{int} dimension of the matrix"},
                  {"howmany=", "{int} number of fourier transforms"},
                  {"niter=", "(int) number of times to repeat the computation"}}
                 );    


    PROFILE_START("Initialize")

    auto dim = args.value<int>("dim", 40);
    auto niter = args.value<int>("niter", 10000);
    auto howmany = args.value<int>("howmany", 4);

    mdarray<std::complex<double>,2> fx({dim*dim, howmany});
    mdarray<std::complex<double>,2> fk({dim*dim, howmany});

    fx.fill(0.);
    fk.fill(0.);

    FourierTransform ft(fx, fk, {dim,dim,1});

    PROFILE_STOP("Initialize")
    

    for(int icounter=0; icounter<niter; ++icounter) {
        PROFILE_START("fft");
        ft.fft(+1);
        PROFILE_STOP("fft");

        PROFILE_START("ifft");
        ft.fft(-1);
        PROFILE_STOP("ifft");
    }

    print_timing(1);
}