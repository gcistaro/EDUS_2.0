#include "Simulation/Simulation.hpp"
#include "core/print_timing.hpp"

int main()
{   
    PROFILE_START("test_simulation")
    auto Rgrid = std::array<int,3>({40,40,1});
    Simulation simulation("/home/gcistaro/NEGF/tb_models/hBN_gap7.25eV_a2.5A", Rgrid);

    simulation.laser.set_Intensity(1.e+05, Wcm2);
    simulation.laser.set_Lambda(800, NanoMeters);
    simulation.laser.set_NumberOfCycles(10);
    simulation.RK_object.set_ResolutionTime(0.001);//0.0003);

    while(simulation.RK_object.get_CurrentTime() < 0.1) {
        simulation.Propagate();
    }
    PROFILE_STOP("test_simulation")
    print_timing(1);
}
