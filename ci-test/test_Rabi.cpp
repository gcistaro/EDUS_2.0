#include "Simulation/Simulation.hpp"
#include "core/print_timing.hpp"
#include "initialize.hpp"

int main()
{   
    initialize();
    PROFILE_START("test_simulation");
    auto Rgrid = std::array<int,3>({1,1,1});
    //Simulation simulation("/home/gcistaro/NEGF/tb_models/hBN_gap7.25eV_a2.5A", Rgrid);
    Simulation simulation("/home/gcistaro/NEGF/tb_models/Rabi2D", Rgrid);

    simulation.laser.set_Intensity(1.e+07, Wcm2);
    //simulation.laser.set_Lambda(3000, NanoMeters);
    simulation.laser.set_Omega(4., ElectronVolt);
    simulation.laser.set_NumberOfCycles(40);
    simulation.RK_object.set_ResolutionTime(simulation.laser.get_Duration()/simulation.laser.get_NumberOfCycles()/1000);

    auto index_FinalTime = simulation.laser.get_Duration()/simulation.RK_object.get_ResolutionTime() + 1;
    for(int it=0; it<index_FinalTime; ++it) {//simulation.RK_object.get_CurrentTime() < 100) {
        simulation.Propagate();
    }
    PROFILE_STOP("test_simulation");
    print_timing(1);
}
