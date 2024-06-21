#include <math.h>
#include "ostream.hpp"

#include "Constants.hpp"
#include "MeshGrid/MeshGrid.hpp"
#include "fftPair/fftPair.hpp"
#include "kGradient/kGradient.hpp"

/* we use the fourier transform of 
     f(i) = std::exp(-a*i*i)
   given by 
     F(w) = std::exp(-pi*pi*w*w/a)   
   the derivative with respect to x is:
   df(x)/dx = -2*a*x*f(x)
      
*/
int main()
{
    /*-----------------------------------------------------------------------------*/
    /*                    this part is a copy of test_fft                          */
    /*           but care, here x is the continuous like variable                  */
    /*-----------------------------------------------------------------------------*/
    size_t Npoints = 1000;
    mdarray<std::complex<double>,2> Array_x({1,Npoints});
    mdarray<std::complex<double>,2> Array_k({1,Npoints});
    std::vector<int> dimensions = {int(Npoints)};
    
    mdarray<double, 1> x({Npoints});
    double i0 = 0;
    double a = Npoints/2.;
    std::ofstream os_x("Array_x.txt");
    for(int i=0; i<Npoints; i++){
        if(i<=Npoints/2){
            x(i) = double(i)/double(Npoints);
        }
        else{
            x(i) = (double(i)-double(Npoints))/double(Npoints);
        }
        Array_x(0,i) = std::exp(-a*std::pow(x(i), 2));
        os_x <<  x(i) << " " <<  Array_x(0,i) << std::endl;
    }
    os_x.close();

    std::cout << "Setting up Fourier Transform..\n";
    FourierTransform fftm(Array_x, Array_k, dimensions);

    std::cout << "Doing Fourier transform...\n";
    fftm.fft(-1);
    std::cout << "DONE!\n";

    /* print array of fourier trasnform */
    std::ofstream os_k("Array_k.txt");
    for(int i=0; i<Npoints; i++){
        os_k << Array_k(0,i) << std::endl;
    }
    os_k.close();
    /*-------------------------------------------------------------------------------------------*/
    /*                        end of the part that is a copy of test_fft                         */

    //getting trivial basis to set up coordinates
    Matrix<double> TrivialMatrix(3,3);
    TrivialMatrix(0,0) = 1.;    TrivialMatrix(0,1) = 0.;    TrivialMatrix(0,2) =  0.;
    TrivialMatrix(1,0) = 0.;    TrivialMatrix(1,1) = 1.;    TrivialMatrix(1,2) =  0.;
    TrivialMatrix(2,0) = 0.;    TrivialMatrix(2,1) = 0.;    TrivialMatrix(2,2) =  1.;

    Basis basis(TrivialMatrix);
    Coordinate::add_Basis(basis, LatticeVectors(R));

    //setting up meshgrid in Fourier space
    std::cout << "setting up meshgrid in Fourier space...\n";
    MeshGrid mg_fft(R, std::array<int,3>({int(Npoints), 1, 1}));

    std::cout << "mg_fft: \n" << mg_fft; 
    std::cout << "done\n";
    //build up gradient object using that meshgrid
    kGradient kgradient;
    kgradient.initialize(mg_fft);

    auto xDerivative_k = Array_k; //to initialize dimensions
    auto xDerivative_x = Array_k;
    Coordinate direction(1,0,0); //direction where we evaluate the gradient
    kgradient.Calculate(xDerivative_k, Array_k, 
                        direction, true);

    std::cout << "Setting up Fourier Transform for derivative..\n";
    FourierTransform fftm_Derivative(xDerivative_x, xDerivative_k, dimensions);

    std::cout << "Doing Fourier transform...\n";
    fftm_Derivative.fft(+1);
    std::cout << "DONE!\n";

    std::cout << "+------------+----------------------+-----------------------+-------------------+\n";
    std::cout << "|    x       | Numerical derivative | Analytical derivative |      Error(%)     |\n";
    std::cout << "+------------+----------------------+-----------------------+-------------------+\n";

    auto&& NumericalSolution = fftm_Derivative.get_Array_x();

    std::ofstream os_NS("Numerical.txt");
    for(int i=0; i<Npoints; i++){
        double xi =x(i);
        //if(i<= Npoints/2.) {
        //    xi = x(i);
        //}else{
        //    xi = x(i)-Npoints;
        //}
        auto AnalyticalSolution = (-2*a*(xi)*Array_x(0,i)).real();
        std::cout << "  |  ";
        std::cout << std::setw(6) << std::setprecision(2) <<  std::scientific << xi;
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
