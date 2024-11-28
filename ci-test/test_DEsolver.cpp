#include <iostream> 
#include <iomanip>
#include <vector>
#include <math.h>
#include <Constants.hpp>
#include "core/print_timing.hpp"
#include "DESolver/DESolver.hpp"
#include "core/cmd_args/cmd_args.hpp"

int main(int argn, char** argv)
{
    cmd_args args(argn, argv, 
                  {{"solver=", "{string} type of solver. supported: RK4, AB4"}}
                 );    
    std::vector<double> Function_;

    auto solver = args.value<std::string>("solver", "RK4");
    auto solver_type = ( solver == "RK4" ? RK4 : AB4 );
    auto threshold   = ( solver == "RK4" ? 1.e-08 : 1.e-04 );

    double InitialTime = 0.;
    double ResolutionTime = 0.001;
    double FinalTime = 10000.;
    static double InitialConstant = 1.;
    static double RateOfIncrease = 2.;

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


    auto DEsolver = DESolver<std::vector<double>>(Function_, InitialCondition, SourceTerm, solver_type);
    
    auto a = RateOfIncrease*InitialConstant*exp(RateOfIncrease*DEsolver.get_CurrentTime() - 0*ResolutionTime);
    auto b = RateOfIncrease*InitialConstant*exp(RateOfIncrease*DEsolver.get_CurrentTime() - 1*ResolutionTime);
    auto c = RateOfIncrease*InitialConstant*exp(RateOfIncrease*DEsolver.get_CurrentTime() - 2*ResolutionTime);
    auto d = RateOfIncrease*InitialConstant*exp(RateOfIncrease*DEsolver.get_CurrentTime() - 3*ResolutionTime);
    DEsolver.set_aux_Function({a}, {b}, {c}, {d});

    DEsolver.set_InitialTime(InitialTime);
    DEsolver.set_ResolutionTime(ResolutionTime);
    


    std::cout << "+---------+------------+--------------------+---------------------+-------------------+\n";
    std::cout << "|  step   |    time    | Numerical solution | Analytical solution |      Error(%)     |\n";
    std::cout << "+---------+------------+--------------------+---------------------+-------------------+\n";

    PROFILE_START("Solve_DE");
    for(double it=0; it<=10000; it++) {
        if(int(it)%100==0){
            auto AnalyticalSolution = InitialConstant*exp(RateOfIncrease*DEsolver.get_CurrentTime());
            auto&& NumericalSolution = DEsolver.get_Function()[0];
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
            auto RelativeError = 100*abs(NumericalSolution-AnalyticalSolution)/abs(AnalyticalSolution);
            std::cout << std::setw(15) << std::setprecision(8) << std::scientific <<  RelativeError;
            std::cout << "  |" << std::endl;
            if( RelativeError > threshold){
                exit(1);
            }
        }
        DEsolver.Propagate();
    }
    PROFILE_STOP("Solve_DE");


    print_timing(1);

    

}
