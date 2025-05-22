#include "Hamiltonian.hpp"
#include "Davidson/davidson.hpp"
#include "core/hdf5/hdf5_tree.hpp"

void BS_Hamiltonian::initialize(std::shared_ptr<Simulation_parameters> ctx__, const int& iq__)
{
    ctx_ = ctx__;
    simulation_ = Simulation(ctx_);

    /* get variables needed */
    auto num_kpoints     = simulation_.DensityMatrix_.get_Operator(Space::k).get_MeshGrid()->get_TotalSize();
    auto num_bands       = simulation_.DensityMatrix_.get_Operator(Space::k).get_nrows();
    auto num_filledbands = ctx_->cfg().filledbands();
    auto num_emptybands  = num_bands - num_filledbands;
    auto kmesh           = simulation_.DensityMatrix_.get_Operator(Space::k).get_MeshGrid();

    
    TwoBodyIndex_.initialize({num_kpoints, num_filledbands, num_emptybands});
    auto Hamiltonian_size = TwoBodyIndex_.get_TotalSize();

    output::title("Problem dimensions");
    output::print("#kpoints           : ", num_kpoints);
    output::print("#bands             : ", num_bands);
    output::print("#filled            : ", num_filledbands);
    output::print("#empty             : ", num_emptybands);
    std::stringstream to_print; 
    to_print << " x " << Hamiltonian_size;
    output::print("Hamiltonian_size   : ", Hamiltonian_size, to_print.str());

    H_BSE_.initialize( Hamiltonian_size, Hamiltonian_size );
    H_BSE_.fill(0.);

    /* IP: diagonal energies */
    auto& Energy = simulation_.Band_energies_;

    #pragma omp parallel for
    for ( int ib=0; ib < Hamiltonian_size; ++ib ) {
        auto nDindex = TwoBodyIndex_.nDindex(ib);
        auto& ik    = nDindex[0];
        auto& ival  = nDindex[1];
        auto& icond = nDindex[2];

        H_BSE_( ib, ib ) += Energy[ik](num_filledbands + icond) - Energy[ik](ival);
    }

    /* off-diagonal energies with Screened potential */
    Operator<std::complex<double>> W_operator; 
    auto Screened_Potential = simulation_.coulomb_.get_ScreenedPotential();
    W_operator.initialize_fft(simulation_.DensityMatrix_, "W");
    std::copy(Screened_Potential.begin(), Screened_Potential.end(), W_operator.get_Operator_R().begin());
    W_operator.go_to_k();
    /* get Matrix of W summed over all R */
    #pragma omp parallel for
    for( int ik1 = 0; ik1 < num_kpoints; ++ik1 ) {
        for( int ik2 = 0; ik2 < num_kpoints; ++ik2 ) {
            /* get variables needed */
            auto q = (*kmesh)[iq__];
            auto k1 = (*kmesh)[ik1];
            auto k2 = (*kmesh)[ik2];
            auto ik2mk1 = kmesh->find( k2-k1 );
            auto ik1pq  = kmesh->find( k1+q ); 
            auto ik2pq  = kmesh->find( k2+q ); 

            auto& Uk1    = Operator<std::complex<double>>::EigenVectors[ik1];
            auto& Uk2    = Operator<std::complex<double>>::EigenVectors[ik2];
            auto& Uk1pq  = Operator<std::complex<double>>::EigenVectors[ik1pq];
            auto& Uk2pq  = Operator<std::complex<double>>::EigenVectors[ik2pq];
            auto& Udk1   = Operator<std::complex<double>>::EigenVectors_dagger[ik1];
            auto& Udk2   = Operator<std::complex<double>>::EigenVectors_dagger[ik2];
            auto& Udk1pq = Operator<std::complex<double>>::EigenVectors_dagger[ik1pq];
            auto& Udk2pq = Operator<std::complex<double>>::EigenVectors_dagger[ik2pq];
            auto detA  = simulation_.jacobian(Coordinate::get_Basis(LatticeVectors(R)).get_M());
            for( int ival1 = 0; ival1 < num_filledbands; ++ival1 ) {
                for( int icond1 = num_filledbands; icond1 < num_bands; ++icond1 ) {
                    auto irow = TwoBodyIndex_.oneDindex(ik1, ival1, icond1-num_filledbands);
                    for( int ival2 = 0; ival2 < num_filledbands; ++ival2 ) {
                        for(int icond2 = num_filledbands; icond2 < num_bands; ++icond2 ) {
                            auto icol = TwoBodyIndex_.oneDindex(ik2, ival2, icond2-num_filledbands);
                            for( int ialpha=0; ialpha < num_bands; ++ialpha ) {
                                for ( int ibeta=0; ibeta < num_bands; ++ibeta ) {
                                    H_BSE_(irow, icol) -= 1./num_kpoints*Udk1pq(icond1, ialpha)*Udk2(ival2, ibeta)*
                                                                          Uk2pq(ialpha, icond2)*Uk1(ibeta, ival1)*
                                                                         W_operator.get_Operator_k()[ik2mk1](ibeta, ialpha);//(ialpha, ibeta);
                                    // == H_BSE_(irow, icol) -= 1./num_kpoints*Udk1pq(ialpha, icond1)*Udk2(ibeta, ival2)*
                                    // ==                                      Uk2pq(icond2, ialpha)*Uk1(ival1, ibeta)*
                                    // ==                                     W_operator.get_Operator_k()[ik2mk1](ialpha, ibeta);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    output::print("Initialize davidson");
    Davidson davidson(&H_BSE_(0,0), Hamiltonian_size, 1, 100);
    output::print("Solve");
    davidson.solve();
/*
    Matrix<std::complex<double>> U;//(Hamiltonian_size, Hamiltonian_size);
    mdarray<double, 1> E;//(std::array<int,1>(Hamiltonian_size));
    H_BSE_.diagonalize(U, E);
    std::ofstream out; 
    out.open("Energy.txt");
    Convert_iterable(E, AuEnergy, ElectronVolt);
    for(int icol=0; icol < Hamiltonian_size; ++icol) {
        out << E(icol) << "     " <<std::endl ;
    }
    out.close();

#ifdef EDUS_HDF5
    HDF5_tree fout("output.h5", hdf5_access_t::read_write);
    fout.create_node(-1);
    fout[-1].create_node("Avc");
    fout[-1]["Avc"].create_node(0);
    
    fout[-1]["Avc"][0].write("local", 
              reinterpret_cast<double*>(&U(0,0)), (U.get_TotalSize() * 2) );
#endif
*/

// == out.open("H_BSE.txt");
// == for(int irow=0; irow < Hamiltonian_size; ++irow) {
// == for(int icol=0; icol < Hamiltonian_size; ++icol) {
// ==     out << H_BSE_(irow, icol).real() << " " << H_BSE_(irow, icol).imag() << "     "; 
// == }
// == out << std::endl; 
// == }
// == 
// == out.close();


}
