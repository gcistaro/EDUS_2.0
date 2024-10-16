#include <iostream> 
#include <iomanip>
#include <vector>
#include <math.h>
#include "RungeKutta/RungeKutta.hpp"
//#include "Adams-Bashforth/Adams-Bashforth.hpp"

//here we solve y'=c*y
//The analytical solution is y=y0*Exp(c*t)

int main()
{
    double InitialTime = 0.;
    double ResolutionTime = 0.001;
    static double InitialConstant = 1.;
    static double RateOfIncrease = 2.;
    std::vector<double> Function_;

    auto InitialCondition = [&](auto& Function){
        Function.resize(1);
        std::fill(Function.begin(), Function.end(), InitialConstant);
    };


    auto SourceTerm = [&](auto& Output, const double InputTime, const auto& InputFunction){
        for(struct {std::vector<double>::iterator Output; std::vector<double>::iterator Input;}
             loop = {Output.begin(), (const_cast<std::vector<double>&>(InputFunction)).begin()};
        (loop.Output!=Output.end()) && (loop.Input != InputFunction.end());  
          ++loop.Output, ++loop.Input)
        {
            (*loop.Output) = RateOfIncrease*(*loop.Input);
        }
    };

    //auto rungekutta = RungeKutta<std::vector<double>>(Function_, InitialCondition, SourceTerm);
    auto rungekutta = RungeKutta<std::vector<double>>(Function_, InitialCondition, SourceTerm);
    rungekutta.set_InitialTime(InitialTime);
    rungekutta.set_ResolutionTime(ResolutionTime);

    std::cout << "+---------+------------+--------------------+---------------------+-------------------+\n";
    std::cout << "|  step   |    time    | Numerical solution | Analytical solution |      Error(%)     |\n";
    std::cout << "+---------+------------+--------------------+---------------------+-------------------+\n";
    for(double it=0; it<=10000; it++){
        if(int(it)%100==0){
            auto AnalyticalSolution = InitialConstant*exp(RateOfIncrease*rungekutta.get_CurrentTime());
            auto&& NumericalSolution = rungekutta.get_Function()[0];
            std::cout << "|";
            std::cout << std::setw(7) << std::fixed << int(it);
            std::cout << "  |  ";
            std::cout << std::setw(6) << std::setprecision(2) <<  std::scientific << rungekutta.get_CurrentTime();
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

        rungekutta.Propagate();
    }

    

}