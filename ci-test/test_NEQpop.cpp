#include "Simulation/Simulation.hpp"
#include "core/print_timing.hpp"
/*
 This test is used to test only the H0 term (no laser). We solve the differential equation:
 i dP(R)/dt = [H,P]
with:
       |     -a             0          |
H(k) = |                               |       (with R fixed)
       |      0             a          |  

          |  0            P0   |
P0(k) =   |                    |
          |  P0*          0    |

that does not commute with the Hamiltonian.
The commutator is: (H.P-P.H)_{01} = -2aP 
Only the off diagonal terms change. The analytical evolution of the 01 term is:

P(k,t) = P0*Exp(i 2*a t)

This is a good proof that the R evolution is working at equilibrium.
*/



int main()
{            
    Simulation simulation("/home/gcistaro/NEGF/tb_models/2B_trivialH", 20.);//;/TBgraphene",40.);//
    
    //modify initialcondition, if not it is 1 for R=0, 0 otherwise
    std::complex<double> a = 0.4+im*0.2;
    std::function<void(BlockMatrix<std::complex<double>, R>&)> InitialConditionToUse = [&](BlockMatrix<std::complex<double>, R>& DM)
    {
        simulation.DensityMatrix.get_Operator_k().fill(0.);
        for(int ik=0; ik < simulation.DensityMatrix.get_Operator_k().get_MeshGrid()->get_TotalSize(); ++ik){
            simulation.DensityMatrix.get_Operator_k()(ik,0,1) = a;
            simulation.DensityMatrix.get_Operator_k()(ik,1,0) = std::conj(a);
        }
        simulation.DensityMatrix.go_to_R();
    };
    //including aliases for sourceterm
    auto& laser = simulation.laser;
    auto& H = simulation.H;
    auto Calculate_TDHamiltonian = [&](const double& time){
        return simulation.Calculate_TDHamiltonian(time);
    };
    #include "Simulation/Functional_SourceTerm.hpp"
    simulation.RK_object.initialize(simulation.DensityMatrix.get_Operator_R(), 
                                    InitialConditionToUse, SourceTerm);

    for(int it=0; it < 10000; ++it){
        simulation.Propagate();
        simulation.DensityMatrix.go_to_k();
        auto DMk = simulation.DensityMatrix.get_Operator_k();
        for(int ik=0; ik < DMk.get_MeshGrid()->get_TotalSize(); ik++){
            //auto& kcart = (*(DMk.get_MeshGrid()))[ik].get("Cartesian");
            //auto& Rcart = (*(simulation.material.H.get_Operator_R().get_MeshGrid()))[1].get("Cartesian");
            //std::cout << "kcart: "<< kcart << "Rcart: " << Rcart;
            //auto dotprod = Rcart[0]*kcart[0]+Rcart[1]*kcart[1]+Rcart[2]*kcart[2];
            auto t = simulation.RK_object.get_CurrentTime();
            std::cout << t << " " << std::setw(20) << std::setprecision(10) << DMk[ik](0,1).real();
            std::cout << std::setw(20) << std::setprecision(10) <<  (a*std::exp(im*2.*simulation.material.H.get_Operator_R()[0](1,1)*t)).real();
            std::cout << std::setw(20) << std::setprecision(10) << DMk[ik](0,1).imag();
            std::cout << std::setw(20) << std::setprecision(10) << (a*std::exp(im*2.*simulation.material.H.get_Operator_R()[0](1,1)*t)).imag();
            std::cout  << std::setw(20) << std::setprecision(10) << std::abs(DMk[ik](0,1)-a*std::exp(im*2.*simulation.material.H.get_Operator_R()[0](1,1)*t))/std::abs(DMk[ik](0,1))*100 << std::endl;
            //if( std::abs( DMk[ik](0,1) - a*std::exp(-im*(1-2*a)*dotprod) ) > 1.e-05){
            //    exit(1);
            //}
        }
    }
}
