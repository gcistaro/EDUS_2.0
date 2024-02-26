#include <math.h>
#include "Constants.hpp"
#include "fftPair/fftPair.hpp"


/* we use the fourier transform of 
     std::exp(-a*i*i)
   given by 
     std::exp(-pi*pi*w*w/a)    
     
*/
int main()
{
    size_t Npoints = 1000;
    mdarray<std::complex<double>,2> Array_x({1,Npoints});
    mdarray<std::complex<double>,2> Array_k({1,Npoints});
    std::vector<int> dimensions = {int(Npoints)};
    
    double i0 = 0;
    double a = 1./double(Npoints);
    for(int i=0; i<Npoints; i++){
        double x;
        if(i<=Npoints/2){
            x = i;
        }
        else{
            x = Npoints-i;
        }
        Array_x(0,i) = std::exp(-a*std::pow(x, 2));
        std::cout << i << " " << x << " " <<  Array_x(0,i) << std::endl;
    }

    std::cout << "Setting up Fourier Transform..\n";
    FourierTransform fftm(Array_x, Array_k, dimensions);

    std::cout << "Doing Fourier transform...\n";
    fftm.fft(-1);
    std::cout << "DONE!\n";

    double DeltaW = 1./double(Npoints);
    double maxW = 1./1.;

    std::cout << "+------------+--------------------+---------------------+-------------------+\n";
    std::cout << "|    freq    | Numerical solution | Analytical solution |      Error(%)     |\n";
    std::cout << "+------------+--------------------+---------------------+-------------------+\n";

    auto&& NumericalSolution = fftm.get_Array_k();

    for(int i=0; i<Npoints; i++){
        double w;
        if(i<=Npoints/2){
            w = double(i)*DeltaW;
        }
        else{
            w = double(Npoints-i)*DeltaW;

        }
        auto AnalyticalSolution = std::sqrt(pi/a)*std::exp(-pi*pi*w*w/a);
        std::cout << "|";
        std::cout << std::setw(7) << std::fixed << int(i);
        std::cout << "  |  ";
        std::cout << std::setw(6) << std::setprecision(2) <<  std::scientific << w;
        std::cout << "  ";
        std::cout << "|";
        std::cout << "  ";
        std::cout << std::setw(16) << std::setprecision(8) << std::scientific <<  NumericalSolution(0,i);
        std::cout << "  ";
        std::cout << "|";
        std::cout << "  ";
        std::cout << std::setw(16) << std::setprecision(8) << std::scientific << AnalyticalSolution;
        std::cout << "   ";
        std::cout << "|";
        std::cout << "  ";
        std::cout << std::setw(15) << std::setprecision(8) << std::scientific <<  100*abs(NumericalSolution(0,i)-AnalyticalSolution)/abs(AnalyticalSolution);
        std::cout << "  |" << std::endl;

        if(AnalyticalSolution > 1.e-07 && 100*abs(NumericalSolution(0,i)-AnalyticalSolution)/abs(AnalyticalSolution)>1.e-04){
            exit(1);
        }
    }

}
