#include "Simulation/Simulation.hpp"
#include "core/print_timing.hpp"
/*
 This test is used to test only the H0 term (no laser). We solve the differential equation:
 i dP(R)/dt = [H,P]
with:
       |     -a(k)                   0          |
H(k) = |                                        |  
       |      0                      a(k)       |  

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
    Simulation simulation("/home/gcistaro/NEGF/tb_models/2B_trivialH", 60.);//;/TBgraphene",40.);//
    
    //modify initialcondition, if not it is 1 for R=0, 0 otherwise
    std::complex<double> a = 0.4+im*0.2;
    std::function<void(BlockMatrix<std::complex<double>, R>&)> InitialConditionToUse = [&](BlockMatrix<std::complex<double>, R>& DM)
    {
        simulation.DensityMatrix.get_Operator_k().fill(0.);
        for(int ik=0; ik < simulation.DensityMatrix.get_Operator_k().get_MeshGrid()->get_TotalSize(); ++ik){
            auto k = (*(simulation.DensityMatrix.get_Operator_k().get_MeshGrid()))[ik].get("LatticeVectors");
            simulation.DensityMatrix.get_Operator_k()(ik,0,1) = std::cos(2.*pi*k[0]);
            simulation.DensityMatrix.get_Operator_k()(ik,1,0) = std::cos(2.*pi*k[0]);
        }
        simulation.DensityMatrix.go_to_R();
    };

    auto DMk0 = simulation.DensityMatrix.get_Operator_k();
    //including aliases for sourceterm
    auto& laser = simulation.laser;
    laser.set_Intensity(0., Wcm2);//1.e+16, Wcm2);
    auto& H = simulation.H;
    auto Calculate_TDHamiltonian = [&](const double& time){
        return simulation.Calculate_TDHamiltonian(time);
    };
    #include "Simulation/Functional_SourceTerm.hpp"
    simulation.RK_object.initialize(simulation.DensityMatrix.get_Operator_R(), 
                                    InitialConditionToUse, SourceTerm);
    simulation.RK_object.set_ResolutionTime(0.001);
    for(int it=0; it < 1000; ++it){
        simulation.Propagate();
        simulation.DensityMatrix.go_to_k();
        auto DMk = simulation.DensityMatrix.get_Operator_k();
        for(int ik=0; ik < DMk.get_MeshGrid()->get_TotalSize(); ik++){
            //auto& kcart = (*(DMk.get_MeshGrid()))[ik].get("Cartesian");
            auto& kfrac = (*(DMk.get_MeshGrid()))[ik].get("LatticeVectors");
            //auto& Rcart = (*(simulation.material.H.get_Operator_R().get_MeshGrid()))[1].get("Cartesian");
            //std::cout << "kcart: "<< kcart << "Rcart: " << Rcart;
            //auto dotprod = Rcart[0]*kcart[0]+Rcart[1]*kcart[1]+Rcart[2]*kcart[2];
            auto t = simulation.RK_object.get_CurrentTime();
            auto Hk = simulation.Band_energies[ik](1);//(0.2+std::cos(2.*pi*kfrac[0]));
            auto Analytical = DMk0[ik](0,1)*std::exp(im*2.*pi*Hk*t);
            auto RelativeError =  std::abs( DMk[ik](0,1) - Analytical)/std::abs(Analytical)*100.;
            std::cout  << std::setw(20) << std::setprecision(10) << it;
            std::cout  << std::setw(20) << std::setprecision(10) << std::abs(DMk[ik](0,1)-Analytical)/std::abs(DMk[ik](0,1))*100 << std::endl;
            if( std::abs(Analytical) > 1.e-07 && 
                RelativeError > 1.){
                exit(1);
            }
        }
    }
}
