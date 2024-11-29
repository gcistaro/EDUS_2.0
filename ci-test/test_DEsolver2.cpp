#include <iostream> 
#include <iomanip>
#include <vector>
#include <math.h>
#include "ostream.hpp"
#include "ConvertUnits.hpp"
#include "Constants.hpp"
#include "DESolver/DESolver.hpp"

//here we solve y'=2*i*H00*y
//to compare with NEQpop2 (same EOM but with numbers!)
//The analytical solution is y=y0*Exp(2*i*H00*t)

int main()
{
    int N = 3;
    double InitialTime = 0.;
    double ResolutionTime = 0.001;
    static double H00 = Convert(2.,ElectronVolt, AuEnergy);
    std::cout << "H00 "<<H00 << std::endl; 
    std::vector<std::complex<double>> Function_;

    auto InitialCondition = [&](auto& Function){
        Function.resize(N);
        for( int ik=0; ik<int(Function.size()); ik++ ) {
            Function[ik] = std::cos(2.*pi*ik/N)/std::sqrt(N);
        }
    };

    auto SourceTerm = [&](auto& Output, const double InputTime, const auto& InputFunction){
        for(struct {std::vector<std::complex<double>>::iterator Output; std::vector<std::complex<double>>::iterator Input;}
             loop = {Output.begin(), (const_cast<std::vector<std::complex<double>>&>(InputFunction)).begin()};
        (loop.Output!=Output.end()) && (loop.Input != InputFunction.end());  
          ++loop.Output, ++loop.Input)
        {
            (*loop.Output) = 2.*im*H00*(*loop.Input);
        }
    };

    std::vector<std::complex<double>> Function0;
    InitialCondition(Function0);

    auto DEsolver = DESolver<std::vector<std::complex<double>>>(Function_, InitialCondition, SourceTerm, RK, 4);
    DEsolver.set_InitialTime(InitialTime);
    DEsolver.set_ResolutionTime(ResolutionTime);

    std::cout << "+---------+------------+--------------------+---------------------+-------------------+\n";
    std::cout << "|  step   |    time    | Numerical solution | Analytical solution |      Error(%)     |\n";
    std::cout << "+---------+------------+--------------------+---------------------+-------------------+\n";
    for(double it=0; it<=1e+4; it++){
        //if(int(it)%100==0){
            for( int ik=0; ik<N; ++ik ) {
                auto AnalyticalSolution = Function0[ik]*std::exp(2.*im*H00*DEsolver.get_CurrentTime());
                auto&& NumericalSolution = DEsolver.get_Function()[ik];
                std::cout << "|";
                std::cout << std::setw(7) << std::fixed << int(it);
                std::cout << "  |  ";
                std::cout << std::setw(6) << std::setprecision(2) <<  std::scientific << DEsolver.get_CurrentTime();
                std::cout << "  ";
                std::cout << "|";
                std::cout << "  ";
                std::cout << std::setw(16) << std::setprecision(8) << std::scientific <<  NumericalSolution;
                std::cout << "  ";
                std::cout << "|";
                std::cout << "  ";
                std::cout << std::setw(16) << std::setprecision(8) << std::scientific << AnalyticalSolution;
                std::cout << "   ";
                std::cout << "|";
                std::cout << "  ";
                std::cout << std::setw(15) << std::setprecision(8) << std::scientific <<  100*abs(NumericalSolution-AnalyticalSolution)/abs(AnalyticalSolution);
                std::cout << "  |" << std::endl;

                if( 100*abs(NumericalSolution-AnalyticalSolution)/abs(AnalyticalSolution) > 1.e-08){
                    exit(1);
                }
            }
        //}    

        DEsolver.Propagate();
    }

    

}