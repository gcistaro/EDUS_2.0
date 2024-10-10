#include "Simulation/Simulation.hpp"
#include "core/print_timing.hpp"
#include "initialize.hpp"
#include "core/projectdir.hpp"

int main()
{   
    initialize();
    PROFILE_START("test_simulation");
    //auto Rgrid = std::array<int,3>({100,100,1});
    //Simulation simulation("/home/gcistaro/NEGF/tb_models/hBN_gap7.25eV_a2.5A", Rgrid);
    std::stringstream json_filename;
    json_filename << ProjectDirectory <<  "/ci-test/inputs/simulation.json";
    Simulation simulation(json_filename.str());

    //simulation.laser.set_Intensity(1.e+10, Wcm2);
    ////simulation.laser.set_Lambda(3000, NanoMeters);
    //simulation.laser.set_Lambda(232.368, NanoMeters);
    ////simulation.laser.set_Omega(10.0, ElectronVolt);
    //simulation.laser.set_NumberOfCycles(1);
    //simulation.RK_object.set_ResolutionTime(0.01);//simulation.laser.get_Duration()/simulation.laser.get_NumberOfCycles()/1000);

    auto index_FinalTime = Convert(150,FemtoSeconds,AuTime)/simulation.RK_object.get_ResolutionTime() +1;
    std::cout << index_FinalTime << std::endl;
    for(int it=0; it<index_FinalTime; ++it) {
        if(it%1000 == 0) {std::cout << "it: " << it << " / " << index_FinalTime<< "  %: " << double(it)/index_FinalTime << std::endl;}
        simulation.Propagate();
    }
    PROFILE_STOP("test_simulation");
    print_timing(1);
}
