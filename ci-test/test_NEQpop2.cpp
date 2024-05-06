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
    auto N = 50;
    Simulation simulation("/home/gcistaro/NEGF/tb_models/2B_trivialH", std::array<int,3>({N,1,1}));//;/TBgraphene",40.);//
    
    std::function<void(BlockMatrix<std::complex<double>, R>&)> InitialConditionToUse = [&](BlockMatrix<std::complex<double>, R>& DM)
    {
        simulation.DensityMatrix.get_Operator_k().fill(0.);
        for(int ik=0; ik < simulation.DensityMatrix.get_Operator_k().get_MeshGrid()->get_TotalSize(); ++ik){
            auto k = (*(simulation.DensityMatrix.get_Operator_k().get_MeshGrid()))[ik].get("LatticeVectors");
            std::cout << N*k[0] << std::endl;
            simulation.DensityMatrix.get_Operator_k()(ik,0,1) = std::cos(2.*pi*k[0])/std::sqrt(N);
            simulation.DensityMatrix.get_Operator_k()(ik,1,0) = std::cos(2.*pi*k[0])/std::sqrt(N);
        }
        simulation.DensityMatrix.go_to_R();
        for(int iR=0; iR < simulation.DensityMatrix.get_Operator_R().get_MeshGrid()->get_TotalSize(); ++iR){
            for(int irow=0; irow<2; ++irow) {
                for(int icol=0; icol<2; ++icol) {
                    if(std::abs(simulation.DensityMatrix.get_Operator_R()[iR](irow, icol)) > 1.e-08 ) {
                        std::cout << iR << " " << irow << " " << icol << " ";
                        std::cout << simulation.DensityMatrix.get_Operator_R()[iR](irow, icol) << std::endl;
                    }
                }
            }
        }
    };


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
    auto DMk0 = simulation.DensityMatrix.get_Operator_k();
    simulation.RK_object.set_ResolutionTime(0.01);

    //std::cout << simulation.DensityMatrix.get_Operator_R();

    auto ST = simulation.DensityMatrix;
    SourceTerm(ST.get_Operator_R(), 0., simulation.DensityMatrix.get_Operator_R());

    std::cout << "SourceTerm\n";
    for(int iR=0; iR < simulation.DensityMatrix.get_Operator_R().get_MeshGrid()->get_TotalSize(); ++iR){
        for(int irow=0; irow<2; ++irow) {
            for(int icol=0; icol<2; ++icol) {
                if(std::abs(ST.get_Operator_R()[iR](irow, icol)) > 1.e-08 ) {
                    std::cout << iR << " " << irow << " " << icol << " ";
                    std::cout << ST.get_Operator_R()[iR](irow, icol) << std::endl;
                }
            }
        }
    }

    std::cout << "Hamiltonian\n";
    for(int iR=0; iR < simulation.material.H.get_Operator_R().get_MeshGrid()->get_TotalSize(); ++iR){
        for(int irow=0; irow<2; ++irow) {
            for(int icol=0; icol<2; ++icol) {
                if(std::abs(simulation.material.H.get_Operator_R()[iR](irow, icol)) > 1.e-08 ) {
                    std::cout << iR << " " << irow << " " << icol << " ";
                    std::cout << simulation.material.H.get_Operator_R()[iR](irow, icol) << std::endl;
                }
            }
        }
    }

    for(int it=0; it < 1000; ++it){
        simulation.Propagate();
        simulation.DensityMatrix.go_to_k();
        auto DMk = simulation.DensityMatrix.get_Operator_k();
        for(int ik=0; ik < DMk.get_MeshGrid()->get_TotalSize(); ik++){
            auto t = simulation.RK_object.get_CurrentTime();
            auto Hk = simulation.Band_energies[ik](1);
            auto Analytical = DMk0[ik](0,1)*std::exp(im*2.*Hk*t);
            auto RelativeError =  std::abs( DMk[ik](0,1) - Analytical)/std::abs(Analytical)*100.;
            std::cout  << std::setw(20) << std::setprecision(10) << it;
            std::cout  << std::setw(40) << std::setprecision(10) << DMk[ik](0,1);
            std::cout  << std::setw(40) << std::setprecision(10) << Analytical;
            std::cout  << std::setw(20) << std::setprecision(10) << RelativeError << std::endl;
            if( std::abs(Analytical) > 1.e-07 && 
                RelativeError > 1.){
                exit(1);
            }
        }
    }
}
