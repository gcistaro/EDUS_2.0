#include "Simulation/Simulation.hpp"
#include "core/print_timing.hpp"

int main()
{            
    Simulation simulation("/home/gcistaro/NEGF/tb_models/hBN_gap7.25eV_a2.5A", 40.);//;/TBgraphene",40.);//
    print_timing(1);
}
