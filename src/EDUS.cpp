

#include <ctime>
#include <iostream>
#include "Simulation/Simulation.hpp"
#include "core/print_timing.hpp"
#include "initialize.hpp"
#include "core/projectdir.hpp"


#define variable(x)  (#x)
#include "core/githash.hpp"



int main(int argc, char *argv[])
{   
    PROFILE_START("EDUS");
    initialize();

    Simulation simulation(argv[1]);

    auto index_FinalTime = Convert(150,FemtoSeconds,AuTime)/simulation.RK_object.get_ResolutionTime() +1;
    std::cout << index_FinalTime << std::endl;
    for(int it=0; it<index_FinalTime; ++it) {
        if(it%1000 == 0) {
            std::cout << "it: " << it << " / " << index_FinalTime<< "  %: "; 
            std::cout << 100*double(it)/index_FinalTime << std::endl;
        }
        simulation.Propagate();
    }
    PROFILE_STOP("test_simulation");
    print_timing(1);
}
