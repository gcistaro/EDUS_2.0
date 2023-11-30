//g++ -std=c++17 -I../header/ test_RK.cpp -o test_RK.x
#include <iostream> 
#include <complex>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <math.h>
#include "RungeKutta.hpp"

//here we solve y'=c*y
//The analytical solution is y=y0*Exp(c*t)

int main()
{
    double InitialTime = 0.;
    double ResolutionTime = 0.01;
    static double InitialConstant = 1.;
    static double RateOfIncrease = 2.;

    auto InitialCondition = [&](auto& Function){
        Function.resize(1);
        std::fill(Function.begin(), Function.end(), InitialConstant);
    };

    auto SourceTerm = [&](auto& Output, const double InputTime, const auto& InputFunction){
        for(auto OutputIterator=Output.begin(), InputFunctionIterator = InputFunction.begin();
        (OutputIterator!=Output.end()) && (InputFunctionIterator != InputFunction.end());  
          ++OutputIterator, ++InputFunctionIterator)
        {
            (*OutputIterator) = RateOfIncrease*(*InputFunctionIterator);
        }
    };

    auto rungekutta = make_RungeKutta<std::vector<double>>(InitialCondition, SourceTerm); 
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
        }    

        rungekutta.Propagate();


    }
    

}