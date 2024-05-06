#include <math.h>
#include "Constants.hpp"
#include "fftPair/fftPair.hpp"


/* we use the fourier transform of 
     std::cos(2*pi*n/N)/std::sqrt(N)
   given by 
     0.5*(delta(1)+delta(-1))    
     
*/
int main()
{
    size_t Npoints = 100;
    mdarray<std::complex<double>,2> Array_x({1,Npoints});
    mdarray<std::complex<double>,2> Array_k({1,Npoints});
    std::vector<int> dimensions = {int(Npoints)};
    
    double i0 = 0;
    double a = 1./double(Npoints);
    for(int i=0; i<Npoints; i++){
        Array_x(0,i) = std::cos(2.*pi*double(i)/double(Npoints))/std::sqrt(double(Npoints));
        std::cout << i <<  Array_x(0,i) << std::endl;
    }

    auto Array_x0 = Array_x;
    std::cout << "Setting up Fourier Transform..\n";
    FourierTransform fftm(Array_x, Array_k, dimensions);

    std::cout << "Doing Fourier transform...\n";
    fftm.fft(-1);
    fftm.fft(+1);

    for(int i=0; i<Npoints; i++){
        std::cout << std::sqrt(Npoints) << " " << std::sqrt(double(Npoints))<<std::endl;
        std::cout << Array_x0(0,i) << " " << Array_x(0,i) << std::endl;///double(Npoints) << std::endl;
    }
    exit(0);

    std::cout << "DONE!\n";


    std::cout << "+------------+--------------------+---------------------+-------------------+\n";
    std::cout << "|    freq    | Numerical solution | Analytical solution |      Error(%)     |\n";
    std::cout << "+------------+--------------------+---------------------+-------------------+\n";

    auto&& NumericalSolution = fftm.get_Array_k();

    for(int i=0; i<Npoints; i++){
        double w;
        if(i<=Npoints/2){
            w = double(i);
        }
        else{
            w = double(Npoints-i);

        }
        double AnalyticalSolution = 0.;
        if (w==1 || w==-1) AnalyticalSolution = 0.5;
        std::cout << "|";
        std::cout << std::setw(7) << std::fixed << int(i);
        std::cout << "  |  ";
        std::cout << std::setw(6) << std::setprecision(2) <<  std::scientific << w;
        std::cout << "  ";
        std::cout << "|";
        std::cout << "  ";
        std::cout << std::setw(16) << std::setprecision(8) << std::scientific <<  NumericalSolution(0,i)*std::sqrt(Npoints);
        std::cout << "  ";
        std::cout << "|";
        std::cout << "  ";
        std::cout << std::setw(16) << std::setprecision(8) << std::scientific << AnalyticalSolution;
        std::cout << "   ";
        std::cout << "|";
        std::cout << "  ";
        std::cout << std::setw(15) << std::setprecision(8) << std::scientific <<  100*abs(NumericalSolution(0,i)-AnalyticalSolution)/abs(AnalyticalSolution);
        std::cout << "  |" << std::endl;

        //if(AnalyticalSolution > 1.e-07 && 100*abs(NumericalSolution(0,i)*std::sqrt(Npoints)-AnalyticalSolution)/abs(AnalyticalSolution)>5.){
        //    exit(1);
        //}
    }

}
