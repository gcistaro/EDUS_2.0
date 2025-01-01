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



int main()
{   
    PROFILE_START("main");
    initialize();
    //------------------------------------------initialize simulation------------------------------------------------------
    auto N = 100;
    std::stringstream ss ;
    ss << ProjectDirectory << "/tb_models/Trivial_Hamiltonian";
    Simulation simulation(ss.str(), std::array<int,3>({N,N,1}));//;/TBgraphene",40.);//

    //------------------------------------------get wanted initial condition------------------------------------------------
    std::function<void(Operator<std::complex<double>>&)> InitialConditionToUse = [&](Operator<std::complex<double>>& DM)
    {
        simulation.DensityMatrix.get_Operator_k().fill(0.);
        for(int ik=0; ik < simulation.DensityMatrix.mpindex.nlocal; ++ik){
            auto ik_local = simulation.DensityMatrix.mpindex.loc1D_to_glob1D(ik);
            auto k_ = (*(simulation.DensityMatrix.get_Operator_k().get_MeshGrid()))[ik_local].get(LatticeVectors(k));
            simulation.DensityMatrix.get_Operator_k()(ik,0,1) = std::exp(2.*pi*im*k_[0])/std::sqrt(N);
            simulation.DensityMatrix.get_Operator_k()(ik,1,0) = std::exp(-2.*pi*im*k_[0])/std::sqrt(N);
        }
        simulation.DensityMatrix.lock_space(Space::k);
    };

    //---------------------------------reinitialize RK with that initial condition--------------------------------------
    auto& setoflaser = simulation.setoflaser;
    Laser laser;
    laser.set_Intensity(1.e+10, Wcm2);
    laser.set_Lambda(800, NanoMeters);
    laser.set_Polarization(Coordinate(1,0,0));
    laser.set_NumberOfCycles(5);
    setoflaser.push_back(laser);
    auto& H = simulation.H;
    auto& kgradient = simulation.kgradient;
    auto& coulomb = simulation.coulomb;
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
    std::cout << "Initializing DEsolver_DM..\n";
    simulation.DEsolver_DM.initialize(simulation.DensityMatrix, 
                                    InitialConditionToUse, SourceTerm, RK,4);
    simulation.DEsolver_DM.set_ResolutionTime(0.1);

    //------------------------------- FINAL CHECK ON PROPAGATOR -----------------------------------------
    std::cout << "Checking correctness of propagator...\n";
    auto DMk0 = simulation.DensityMatrix.get_Operator_k();
    for(int it=0; it <= 20; ++it){
        simulation.DensityMatrix.go_to_k();
        auto& DMk = simulation.DensityMatrix.get_Operator_k();
        for(int ik_loc=0; ik_loc < simulation.DensityMatrix.mpindex.nlocal; ++ik_loc){
            auto ik_glob = simulation.DensityMatrix.mpindex.loc1D_to_glob1D(ik_loc);
            auto k_ = (*(simulation.DensityMatrix.get_Operator_k().get_MeshGrid()))[ik_glob].get(LatticeVectors(k));
            auto t = simulation.DEsolver_DM.get_CurrentTime();

            auto Analytical = std::exp(2.*pi*im*(k_[0]+laser.VectorPotential(t).get(LatticeVectors(k))[0]))/std::sqrt(N);
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
