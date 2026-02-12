#include <cassert>
#include "initialize.hpp"
#include "Simulation/Simulation.hpp"
#include "core/print_timing.hpp"
#include "core/projectdir.hpp"

/*
 This test is used to test only the adiabatic term.
 We solve:
     dP(k,t)/dt = E(t).Nabla_k P(k,t)
 The result is:
    P(k,t) = P(k-A(t)+A(t0), t0)

 we choose as density matrix: 
    P(k,t0) = [      0             e^{2\pi i k_0} ]
              [e^{-2\pi i k_0}           0        ]
*/



int main(int argc, char *argv[])
{   
    PROFILE_START("main");
    initialize();

    auto ctx = std::make_shared<Simulation_parameters>();
    ctx->import(std::string(argv[1]));
    Simulation simulation(ctx);
  
    /* get aliases for restarting simulation with the new initial conditions */
    auto& DensityMatrix_ = simulation.DensityMatrix_;
    auto& setoflaser_ = simulation.setoflaser_;
    auto& DEsolver_DM_ = simulation.DEsolver_DM_;
    auto& H_ = simulation.H_;
    auto& kgradient_ = simulation.kgradient_;
    auto& coulomb_ = simulation.coulomb_;
    auto& SpaceOfPropagation_Gradient_ = simulation.SpaceOfPropagation_Gradient_;
    auto& ctx_ = simulation.ctx_;
    auto Apply_Peierls_phase =
        [&](auto&&... args) {
            return simulation.Apply_Peierls_phase(
                std::forward<decltype(args)>(args)...);
        };
    auto& SpaceOfPropagation_ = simulation.SpaceOfPropagation_;

    /* function for the new initial condition to replace the one in SBE */
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

    /* reinitialize DE_Solver with the new initial conditions */
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

    /* propagate numerically and compare with analytical solution */
    std::cout << "Checking correctness of propagator...\n";
    Operator<std::complex<double>> DM0;
    DM0.initialize_fft(simulation.DensityMatrix_);
    
    for(int it=0; it <= setoflaser_[0].get_Duration()/DEsolver_DM_.get_ResolutionTime(); ++it){
        DensityMatrix_.go_to_k();
        std::stringstream ss; ss << "dm_"<<it<<".txt";
        std::ofstream os(ss.str());
        auto& DMk = DensityMatrix_.get_Operator_k();
        InitialConditionToUse(DM0);

        for(int ik_loc=0; ik_loc < DensityMatrix_.mpindex.nlocal; ++ik_loc){
            auto ik_glob = DensityMatrix_.mpindex.loc1D_to_glob1D(ik_loc);
            auto k_ = (*(DensityMatrix_.get_Operator_k().get_MeshGrid()))[ik_glob].get(LatticeVectors(k));
            auto t = DEsolver_DM_.get_CurrentTime();

            simulation.Apply_Peierls_phase(DM0, t, -1);
            DM0.go_to_k();
            auto Analytical = DM0.get_Operator_k()(ik_loc,0,1);//std::exp(2.*pi*im*(k_[0]-setoflaser_.VectorPotential(t).get(LatticeVectors(k))[0]));
            
            auto RelativeError =  std::abs( DMk[ik_loc](0,1) - Analytical)/std::abs(Analytical)*100.;
            //os  << std::setw(5) << std::setprecision(10) << k_[0];
            //os  << std::setw(20) << std::setprecision(10) << setoflaser_.VectorPotential(t).get(LatticeVectors(k))[0];
            //os  << std::setw(40) << std::setprecision(10) << DMk[ik_loc](0,1);
            //os  << std::setw(40) << std::setprecision(10) << Analytical;
            //os << std::endl;
            //std::cout  << std::setw(20) << std::setprecision(10) << RelativeError << std::endl;
            if( std::abs(Analytical) > 1.e-07 && 
                RelativeError > 1.e-03){
                exit(1);
            }
        }
        os.close();
        simulation.do_onestep();
    }

    PROFILE_STOP("main");

    print_timing(1);
}
