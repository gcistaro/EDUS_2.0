#include <math.h>
#include "../Constants.hpp"
#include "fftPair.hpp"


/* we use the fourier transform of 
     std::exp(-a*i*i)
   given by 
     std::exp(-pi*pi*w*w/a)    
     
*/
int main()
{
    int Npoints = 1000;
    mdarray<std::complex<double>,1> Input({Npoints});
    
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
        Input(i) = std::exp(-a*std::pow(x, 2));
        std::cout << i << " " << x << " " <<  Input(i) << std::endl;
    }

    fftPair fftm(Input);

    fftm.fft(-1);

    double DeltaW = 1./double(Npoints);
    double maxW = 1./1.;

    std::cout << "+------------+--------------------+---------------------+-------------------+\n";
    std::cout << "|    freq    | Numerical solution | Analytical solution |      Error(%)     |\n";
    std::cout << "+------------+--------------------+---------------------+-------------------+\n";

    auto&& NumericalSolution = fftm.get_OutputArray();

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
        std::cout << std::setw(16) << std::setprecision(8) << std::scientific <<  NumericalSolution(i);
        std::cout << "  ";
        std::cout << "|";
        std::cout << "  ";
        std::cout << std::setw(16) << std::setprecision(8) << std::scientific << AnalyticalSolution;
        std::cout << "   ";
        std::cout << "|";
        std::cout << "  ";
        std::cout << std::setw(15) << std::setprecision(8) << std::scientific <<  100*abs(NumericalSolution(i)-AnalyticalSolution)/abs(AnalyticalSolution);
        std::cout << "  |" << std::endl;
    }

}
