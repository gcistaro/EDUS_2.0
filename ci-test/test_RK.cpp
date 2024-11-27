#include <iostream> 
#include <iomanip>
#include <vector>
#include <math.h>
#include "DESolver/DESolver.hpp"
#include <chrono>

//here we solve y'=c*y
//The analytical solution is y=y0*Exp(c*t)

int main()
{
    double InitialTime = 0.;
    double ResolutionTime = 0.00001;
    double FinalTime = 1.;
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

    /*
    auto rungekutta = RungeKutta<std::vector<double>>(Function_, InitialCondition, SourceTerm);
    //auto rungekutta = DESolver<std::vector<double>>(Function_, InitialCondition, SourceTerm);             // how do i specify that it is the derived class?
    //RungeKutta<decltype(Function_)> rungekutta(Function_, InitialCondition, SourceTerm);                  // here i define the object rungekutta as a RungeKutta derived class object
    rungekutta.set_InitialTime(InitialTime);
    rungekutta.set_ResolutionTime(ResolutionTime);
    */

    /*
    auto adamsbashforth = AdamsBashforth<std::vector<double>>(Function_, InitialCondition, SourceTerm);
    adamsbashforth.set_InitialTime(InitialTime);
    adamsbashforth.set_ResolutionTime(ResolutionTime);

    */
    

    // uncomment the following lines if you include AdamsBashforth and not RungeKutta

    auto PropagatedFunction = DESolver<std::vector<double>>(Function_, InitialCondition, SourceTerm, DESolver<decltype(Function_)>::AB4);
    
    auto a = RateOfIncrease*InitialConstant*exp(RateOfIncrease*PropagatedFunction.get_CurrentTime() - 0*ResolutionTime);
    auto b = RateOfIncrease*InitialConstant*exp(RateOfIncrease*PropagatedFunction.get_CurrentTime() - 1*ResolutionTime);
    auto c = RateOfIncrease*InitialConstant*exp(RateOfIncrease*PropagatedFunction.get_CurrentTime() - 2*ResolutionTime);
    auto d = RateOfIncrease*InitialConstant*exp(RateOfIncrease*PropagatedFunction.get_CurrentTime() - 3*ResolutionTime);
    PropagatedFunction.setFnsAB4({a}, {b}, {c}, {d});

   PropagatedFunction.set_InitialTime(InitialTime);
   PropagatedFunction.set_ResolutionTime(ResolutionTime);
    

    auto start = std::chrono::high_resolution_clock::now();

    std::cout << "+---------+------------+--------------------+---------------------+-------------------+\n";
    std::cout << "|  step   |    time    | Numerical solution | Analytical solution |      Error(%)     |\n";
    std::cout << "+---------+------------+--------------------+---------------------+-------------------+\n";
    //for(double it=0; it<=10000; it++){
    int it = 0;
    while (PropagatedFunction.get_CurrentTime() <= FinalTime){
        if(int(it)%100==0){
            auto AnalyticalSolution = InitialConstant*exp(RateOfIncrease*PropagatedFunction.get_CurrentTime());
            auto&& NumericalSolution = PropagatedFunction.get_Function()[0];
            std::cout << "|";
            std::cout << std::setw(7) << std::fixed << int(it);
            std::cout << "  |  ";
            std::cout << std::setw(6) << std::setprecision(2) <<  std::scientific << PropagatedFunction.get_CurrentTime();
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

        PropagatedFunction.Propagate();
        it++;
    }

    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> duration = end - start;
    std::cout << "Time taken: " << duration.count() << " seconds" << std::endl;

    

}