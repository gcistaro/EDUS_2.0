#include "Simulation/Simulation.hpp"
#include "core/print_timing.hpp"
/*
 This test is used to test only the R term. We solve the differential equation:
 i dP(R)/dt = (E(t).R) P(R)
 The solution is analytical:
 P(R,t) = Exp(-i R.\int_{t_0}^t E(t')dt')P(R,t_0)
 or with the vector potential:
  P(R,t) = Exp(+i R.A(t))P(R,t_0)
to simplify things we set 
P(R,t0)=1 for every R.
*/



int main()
{            
    Simulation simulation("/home/gcistaro/NEGF/tb_models/Trivial_Hamiltonian", 20.);//;/TBgraphene",40.);//
    
    //modify initialcondition, if not it is 1 for R=0, 0 otherwise
    std::function<void(BlockMatrix<std::complex<double>, R>&)> InitialConditionToUse = [&](BlockMatrix<std::complex<double>, R>& DM)
    {
        for(int iR=0; iR < simulation.DensityMatrix.get_Operator_R().get_MeshGrid()->get_TotalSize(); iR++){
            simulation.DensityMatrix.get_Operator_R()[iR](0,0) = 1.;
        }
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

    auto DMR0 = simulation.DensityMatrix.get_Operator_R();


    //std::ofstream os("Las.txt");
    //for(int it=0; it < 20000; ++it){
    //    auto las = simulation.laser(it*0.1).get("Cartesian");
    //    auto A = simulation.laser.VectorPotential(it*0.1).get("Cartesian");
    //    os << it*0.1 << " " << las[0] << " "<< A[0] << std::endl;
    //}
    //os.close();
    for(int it=0; it < 1000; ++it){
        simulation.Propagate();
        auto DMR = simulation.DensityMatrix.get_Operator_R();
        auto A = simulation.laser.VectorPotential(simulation.RK_object.get_CurrentTime()).get("Cartesian");
        for(int iR=0; iR < DMR.get_MeshGrid()->get_TotalSize(); iR++){
            auto& Rcart = (*(DMR.get_MeshGrid()))[iR].get("Cartesian");
            auto dotprod = Rcart[0]*A[0]+Rcart[1]*A[1]+Rcart[2]*A[2];
            std::cout << simulation.RK_object.get_CurrentTime() << " " << std::abs( DMR[iR](0,0) - std::exp(im*dotprod) ) << std::endl;
            if( std::abs( DMR[iR](0,0) - std::exp(im*dotprod) ) > 1.e-05){
                exit(1);
            }
        }
    }

}
