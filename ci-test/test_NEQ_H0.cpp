#include <cassert>
#include "initialize.hpp"
#include "Simulation/Simulation.hpp"
#include "core/print_timing.hpp"
#include "core/projectdir.hpp"
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
    initialize();
    //------------------------------------------initialize simulation------------------------------------------------------
    auto N = 100;
    std::stringstream ss ;
    ss << ProjectDirectory << "/tb_models/2B_trivialH";
    Simulation simulation(ss.str(), std::array<int,3>({N,N,1}));//;/TBgraphene",40.);//

    //------------------------------------------get wanted initial condition------------------------------------------------
    std::function<void(Operator<std::complex<double>>&)> InitialConditionToUse = [&](Operator<std::complex<double>>& DM)
    {
        simulation.DensityMatrix.get_Operator_k().fill(0.);
        for(int ik=0; ik < simulation.DensityMatrix.mpindex.nlocal; ++ik){
            auto ik_local = simulation.DensityMatrix.mpindex.loc1D_to_glob1D(ik);
            auto k_ = (*(simulation.DensityMatrix.get_Operator_k().get_MeshGrid()))[ik_local].get(LatticeVectors(k));
            simulation.DensityMatrix.get_Operator_k()(ik,0,1) = std::cos(2.*pi*k_[0])/std::sqrt(N);
            simulation.DensityMatrix.get_Operator_k()(ik,1,0) = std::cos(2.*pi*k_[0])/std::sqrt(N);
        }
        simulation.DensityMatrix.lock_space(Space::k);
    };

    //---------------------------------reinitialize RK with that initial condition--------------------------------------
    auto& setoflaser = simulation.setoflaser;
    Laser laser;
    laser.set_Intensity(0, Wcm2);
    laser.set_Lambda(5, NanoMeters);
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
    std::cout << "Initializing RK_object..\n";
    simulation.DEsolver.initialize(simulation.DensityMatrix, 
                                    InitialConditionToUse, SourceTerm, RK, 4);
    simulation.DEsolver.set_ResolutionTime(0.1);
    //---------------------------------check correctness of Source term--------------------------------------
    //Operator<std::complex<double>> ST;
    //ST.initialize_fft(*(simulation.DensityMatrix.get_Operator_R().get_MeshGrid()), 
    //                                    simulation.DensityMatrix.get_Operator_R().get_nrows());
    //std::cout << "ST: " << &ST << std::endl; 
    //ST.go_to_k();
    //SourceTerm(ST, 0., simulation.DensityMatrix);
    //for( int ik=0; ik<ST.get_Operator_k().get_nblocks(); ++ik ) {
    //    auto& H11 = simulation.material.H.get_Operator_k()[ik](1,1);
    //    auto& P01 = simulation.DensityMatrix.get_Operator_k()[ik](0,1);
    //    auto& ST_k = ST.get_Operator_k();
    //    if(std::abs(ST.get_Operator_k()[ik](0,0)) > 1.e-10 || std::abs(ST.get_Operator_k()[ik](1,1)) > 1.e-10) {
    //        std::cout << "ST(0,0) " << ST_k(ik,0,0) <<" or ST(1,1) "<< ST_k(ik,1,1) << " is diffferent than the analytical one!\n";
    //        exit(1);
    //    }
    //    if(std::abs(ST.get_Operator_k()[ik](0,1) - (im*2.*P01*H11))>1.e-10) {
    //        std::cout << "ST(0,1) " << ST_k(ik,0,1) << " is diffferent than the analytical one "<< im*2.*P01*H11<<"!\n";
    //        exit(1);
    //    }
    //}
    //std::cout << "SourceTerm is correctly calculated!!\n";


    //------------------------------- FINAL CHECK ON PROPAGATOR -----------------------------------------
    std::cout << "Checking correctness of propagator...\n";

    auto DMk0 = simulation.DensityMatrix.get_Operator_k();
    //for(int it=0; it <= 20; ++it){
    //    simulation.DensityMatrix.go_to_k();
    //    auto& DMk = simulation.DensityMatrix.get_Operator_k();
    //    for(int ik=0; ik < simulation.DensityMatrix.mpindex.nlocal; ++ik){
    //        auto ik_local = simulation.DensityMatrix.mpindex.loc1D_to_glob1D(ik);
    //        auto k_ = (*(simulation.DensityMatrix.get_Operator_k().get_MeshGrid()))[ik_local].get(LatticeVectors(k));
    //        auto t = simulation.RK_object.get_CurrentTime();
    //        auto Hk = simulation.Band_energies[ik_local](1);
    //        auto Analytical = DMk0[ik_local](0,1)*std::exp(im*2.*Hk*t);
    //        auto RelativeError =  std::abs( DMk[ik](0,1) - Analytical)/std::abs(Analytical)*100.;
//
    //        std::cout  << std::setw(20) << std::setprecision(10) << it;
    //        std::cout  << std::setw(40) << std::setprecision(10) << DMk[ik](0,1);
    //        std::cout  << std::setw(40) << std::setprecision(10) << Analytical;
    //        std::cout  << std::setw(20) << std::setprecision(10) << RelativeError << std::endl;
    //        if( std::abs(Analytical) > 1.e-07 && 
    //            RelativeError > 10.){
    //            exit(1);
    //        }
    //    }
    //    simulation.Propagate();
    //}


    for(int it=0; it <= 20; ++it){
        simulation.DensityMatrix.go_to_k();
        auto DMk = simulation.DensityMatrix.get_Operator_k();
        for(int ik=0; ik < DMk.get_nblocks(); ik++){
            auto t = simulation.DEsolver.get_CurrentTime();
            auto Hk = simulation.Band_energies[ik](1);
            auto Analytical = DMk0[ik](0,1)*std::exp(im*2.*Hk*t);
            auto RelativeError =  std::abs( DMk[ik](0,1) - Analytical)/std::abs(Analytical)*100.;
            std::cout  << std::setw(20) << std::setprecision(10) << it;
            std::cout  << std::setw(20) << std::setprecision(10) << ik;
            std::cout  << std::setw(40) << std::setprecision(10) << DMk[ik](0,1);
            std::cout  << std::setw(40) << std::setprecision(10) << Analytical;
            std::cout  << std::setw(20) << std::setprecision(10) << RelativeError << std::endl;
            if( std::abs(Analytical) > 1.e-07 && 
                RelativeError > 1.e-06){
                exit(1);
            }
        }
        simulation.do_onestep();
    }
    print_timing(1);
}
