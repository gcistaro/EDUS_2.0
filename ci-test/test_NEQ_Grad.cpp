#include <cassert>
#include "initialize.hpp"
#include "Simulation/Simulation.hpp"
#include "core/print_timing.hpp"
#include "core/projectdir.hpp"

/*
 This test is used to test only the adiabatic term.
 We solve:
     idP(k,t)/dt = E(t).Nabla_k P(k,t)
 The result is:
    P(k,t) = P(k+A(t)-A(t0), t0)
*/



int main(int argc, char *argv[])
{   
    PROFILE_START("main");
    initialize();
    //------------------------------------------initialize simulation------------------------------------------------------
    /*
    auto N = 100;
    std::stringstream ss ;
    ss << ProjectDirectory << "/tb_models/Trivial_Hamiltonian";
    */
    auto ctx = std::make_shared<Simulation_parameters>();
    ctx->import(std::string(argv[1]));
    Simulation simulation(ctx);
  
    auto& DensityMatrix_ = simulation.DensityMatrix_;
    auto& setoflaser_ = simulation.setoflaser_;
    auto& DEsolver_DM_ = simulation.DEsolver_DM_;
    auto& H_ = simulation.H_;
    auto& kgradient_ = simulation.kgradient_;
    auto& coulomb_ = simulation.coulomb_;


    //------------------------------------------get wanted initial condition------------------------------------------------
    std::function<void(Operator<std::complex<double>>&)> InitialConditionToUse = [&](Operator<std::complex<double>>& DM)
    {
        DensityMatrix_.get_Operator_k().fill(0.);
        for(int ik=0; ik < DensityMatrix_.mpindex.nlocal; ++ik){
            auto ik_local = DensityMatrix_.mpindex.loc1D_to_glob1D(ik);
            auto k_ = (*(DensityMatrix_.get_Operator_k().get_MeshGrid()))[ik_local].get(LatticeVectors(k));
            DensityMatrix_.get_Operator_k()(ik,0,1) = std::exp(2.*pi*im*k_[0]);
            DensityMatrix_.get_Operator_k()(ik,1,0) = std::exp(-2.*pi*im*k_[0]);
        }
        DensityMatrix_.lock_space(Space::k);
    };

    //---------------------------------reinitialize RK with that initial condition--------------------------------------
    auto Calculate_TDHamiltonian = [&](const double& time, const bool& erase){
        return simulation.Calculate_TDHamiltonian(time, erase);
    };
    auto get_it = [&](const double& time){
        return simulation.get_it(time);
    };
    auto PrintObservables = [&](const double& time){
        return simulation.PrintObservables(time);
    };

    #include "Simulation/Functional_SourceTerm.hpp"
    std::cout << "Initializing DEsolver_DM_..\n";
    DEsolver_DM_.initialize(DensityMatrix_, 
                                    InitialConditionToUse, SourceTerm, RK,4);
    DEsolver_DM_.set_ResolutionTime(0.1);

    //------------------------------- FINAL CHECK ON PROPAGATOR -----------------------------------------
    std::cout << "Checking correctness of propagator...\n";
    auto DMk0 = DensityMatrix_.get_Operator_k();
    for(int it=0; it <= 20; ++it){
        DensityMatrix_.go_to_k();
        auto& DMk = DensityMatrix_.get_Operator_k();
        for(int ik_loc=0; ik_loc < DensityMatrix_.mpindex.nlocal; ++ik_loc){
            auto ik_glob = DensityMatrix_.mpindex.loc1D_to_glob1D(ik_loc);
            auto k_ = (*(DensityMatrix_.get_Operator_k().get_MeshGrid()))[ik_glob].get(LatticeVectors(k));
            auto t = DEsolver_DM_.get_CurrentTime();

            auto Analytical = std::exp(2.*pi*im*(k_[0]+setoflaser_.VectorPotential(t).get(LatticeVectors(k))[0]));
            auto RelativeError =  std::abs( DMk[ik_loc](0,1) - Analytical)/std::abs(Analytical)*100.;
            std::cout  << std::setw(20) << std::setprecision(10) << it;
            std::cout  << std::setw(40) << std::setprecision(10) << DMk[ik_loc](0,1);
            std::cout  << std::setw(40) << std::setprecision(10) << Analytical;
            std::cout  << std::setw(20) << std::setprecision(10) << RelativeError << std::endl;
            if( std::abs(Analytical) > 1.e-07 && 
                RelativeError > 1.e-03){
                exit(1);
            }
        }
        simulation.do_onestep();
    }

    PROFILE_STOP("main");

    print_timing(1);
}
