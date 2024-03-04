#include "Simulation/Simulation.hpp"
#include "core/print_timing.hpp"

int main()
{            
    Simulation simulation("/home/gcistaro/NEGF/tb_models/non_dummy_gap", 20.);
    print_timing(1);
}
