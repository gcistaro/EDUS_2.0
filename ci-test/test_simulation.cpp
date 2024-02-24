#include "Simulation/Simulation.hpp"

int main()
{            
    //Filling DM at t0 with electron for bands below Fermi energy
    auto InitialCondition = [&](Operator<std::complex<double>>& DM){
        std::fill(DM.get_Operator_R().begin(), DM.get_Operator_R().end(), 0.);
    };

    //meaning i dP/dT = [P, H] with H=H0+E.r
    auto SourceTerm = [&](Operator<std::complex<double>>& DM_Output, const double& time, const Operator<std::complex<double>>& DM_Input, const double& time1){
        std::fill(DM_Output.get_Operator_R().begin(), DM_Output.get_Operator_R().end(), 0.);
    };

    Simulation simulation;
}