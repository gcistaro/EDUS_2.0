#include <math.h>
#include "Constants.hpp"
#include "ostream.hpp"

#include "fftPair/fftPair.hpp"


/* we use the fourier transform of 
     f(i) = std::exp(-a*i*i)
   given by 
     F(w) = std::exp(-pi*pi*w*w/a)    
     
*/
int main()
{
    int Npoints = 1000;
    mdarray<std::complex<double>,2> Array_x({1,Npoints});
    mdarray<std::complex<double>,2> Array_k({1,Npoints});
    std::vector<int> dimensions = {int(Npoints)};
    
    double a = 1./double(Npoints);
    std::ofstream os_x("Array_x.txt");
    for(int i=0; i<Npoints; i++){
        double x;
        if(i<=Npoints/2){
            x = i;
        }
        else{
            x = Npoints-i;
        }
        Array_x(0,i) = std::exp(-a*std::pow(x, 2));
        os_x << x << " " <<  Array_x(0,i) << std::endl;
    }
    os_x.close();

    std::cout << "Setting up Fourier Transform..\n";
    FourierTransform fftm(Array_x, Array_k, dimensions);

    std::cout << "Doing Fourier transform...\n";
    fftm.fft(-1);
    std::cout << "DONE!\n";
    std::ofstream os_k("Array_k.txt");
    for(int i=0; i<Npoints; i++){
        os_k << Array_k(0,i) << std::endl;
    }
    os_k.close();

    double DeltaW = 1./double(Npoints);

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
        auto AnalyticalSolution = std::sqrt(pi/a)*std::exp(-pi*pi*w*w/a)/std::sqrt(Npoints);
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

        if(AnalyticalSolution > 1.e-07 && 100*abs(NumericalSolution(0,i)-AnalyticalSolution)/abs(AnalyticalSolution)>5.){
            exit(1);
        }
    }

}
