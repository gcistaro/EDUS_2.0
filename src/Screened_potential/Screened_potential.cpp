#include "Screened_potential.hpp"
#include "omp.h"

void ScreenedPotential::initialize()
{
    /* Initialize all the dimensions */
    W_             .initialize_fft(DensityMatrix_);
    V_             .initialize_fft(DensityMatrix_);
    X_             .initialize_fft(DensityMatrix_);
    Epsilon_       .initialize_fft(DensityMatrix_);
    InverseEpsilon_.initialize_fft(DensityMatrix_);
    nodename_ = "-1";
}

void ScreenedPotential::ResponseFunction()
{
    /* Lindhard formula in point like approximation to get X_{ij}(q) */

    X_.go_to_k(false);
    auto& Xk = X_.get_Operator(Space::k);
    Xk.fill(0.);

    auto num_bands       =  Xk.get_nrows();
    auto num_filledbands =  ctx_->cfg().filledbands();
    auto num_emptybands  =  num_bands - num_filledbands;
    auto kmesh           =  Xk.get_MeshGrid();

    #pragma omp parallel for
    for( int iq = 0; iq < kmesh->get_TotalSize(); iq++ ) {
        output::print("Calculation for iq = ", iq);
        auto q = (*kmesh)[iq];
        for ( int ik = 0; ik < kmesh->get_TotalSize(); ++ik ) {
            auto k = (*kmesh)[ik];
            auto ikp = kmesh->find( k+q );//index of k+q

            auto& Uk   = Operator<std::complex<double>>::EigenVectors[ik];
            auto& Udkp = Operator<std::complex<double>>::EigenVectors_dagger[ikp]; 
            auto& Ukp  = Operator<std::complex<double>>::EigenVectors[ikp]; 
            auto& Udk  = Operator<std::complex<double>>::EigenVectors_dagger[ik];

            for( int ialpha=0; ialpha < num_bands; ++ialpha ) {
                for ( int ibeta=0; ibeta < num_bands; ++ibeta ) {
                    for( int icond = num_filledbands; icond < num_bands; ++icond ) {
                        for ( int ival = 0; ival < num_filledbands; ++ival ) {
                            //sum over k, v, c 
                            Xk(iq, ialpha, ibeta) += 
                                (Uk(ialpha, ival)*Ukp(ibeta, icond)*Udkp(icond, ialpha)*Udk(ival, ibeta)) /
                                                        (kmesh->get_TotalSize()*(Band_energies_[ik](ival)-Band_energies_[ikp](icond)));
                        }
                    }
                }//ibeta
            }//ialpha
        }//ik 
    }//iq

    X_.go_to_R();

    X_.print_Rdecay("X", material_.rwann_);

#ifdef EDUS_HDF5
    std::string name = "output.h5";
    if (!file_exists(name)) {
        HDF5_tree(name, hdf5_access_t::truncate);
    }
    HDF5_tree fout(name, hdf5_access_t::read_write);
    mpi::Communicator::world().barrier();
    fout.create_node(nodename_);
    X_.get_Operator_k().write_h5("output.h5",nodename_,"Xk");
    Operator<std::complex<double>>::EigenVectors.write_h5("output.h5", nodename_, "Uk");
#endif

}

void ScreenedPotential::BareCoulomb()
{
    /* compute V_{ij}(R) = 1./(r_i+R-r_j) */
    V_.go_to_R(false);

    auto& VR = V_.get_Operator(Space::R);
    VR.fill(0.);

    auto Rgamma_centered = get_GammaCentered_grid(*VR.get_MeshGrid());

    for( int iR_loc=0; iR_loc < VR.get_nblocks(); ++iR_loc ) {
        auto iR_glob = iR_loc;//TODO change for mpi
        auto& Rvec = Rgamma_centered[iR_glob];
        for( int irow = 0; irow < VR.get_nrows(); ++irow ) {
            for(int icol = 0; icol < VR.get_ncols(); ++icol) {
                auto ratom = material_.rwann_[irow] + Rvec - material_.rwann_[icol];
                if( ratom.norm() < 1.e-05 ) {
                    VR(iR_loc, irow, icol) = 0.;
                }
                else {
                    VR(iR_loc, irow, icol) = 1./ratom.norm();
                }

            }

        }
    }

    V_.print_Rdecay("V", material_.rwann_);

#ifdef EDUS_HDF5
    V_.go_to_k(true);
    HDF5_tree fout("output.h5", hdf5_access_t::read_write);
    mpi::Communicator::world().barrier();
    V_.get_Operator_k().write_h5("output.h5",nodename_,"Vk");
    V_.go_to_R(false);
#endif
}

void ScreenedPotential::Epsilon()
{
    /* Implementation of Epsilon(r,r') = \delta(r-r') -\int d^3 r'' V(r,r'')X(r'',r') */

    /* -V_{\alpha\beta}X_{\alpha\beta} + 1 */
    Epsilon_.go_to_R(false);

    auto& EpsR = Epsilon_.get_Operator(Space::R);
    auto& VR   = V_      .get_Operator(Space::R);
    auto& XR   = X_      .get_Operator(Space::R);
    auto kmesh           =  XR.get_MeshGrid();


    #pragma omp parallel for
    for( int iR = 0; iR<EpsR.get_nblocks(); ++iR ) {
        for( int irow = 0; irow < EpsR.get_nrows(); ++irow ) {
            for( int icol = 0; icol < EpsR.get_ncols(); ++icol ) {
                EpsR(iR, irow, icol) = -VR(iR, irow, icol)*XR(iR, irow, icol) + 1.;
            }
        }
    }

    Epsilon_.print_Rdecay("Epsilon", material_.rwann_);

#ifdef EDUS_HDF5
    Epsilon_.go_to_k(true);
    HDF5_tree fout("output.h5", hdf5_access_t::read_write);
    mpi::Communicator::world().barrier();
    Epsilon_.get_Operator_k().write_h5("output.h5",nodename_,"Epsilonk");
    Epsilon_.go_to_R(false);
#endif
}

void ScreenedPotential::Calculate()
{
    ResponseFunction();
    BareCoulomb();
    Epsilon();

    /* get inverse of epsilon */
    InverseEpsilon_.go_to_R(false);

    auto& invEpsR = InverseEpsilon_.get_Operator(Space::R);
    auto& EpsR    = Epsilon_       .get_Operator(Space::R);

    /* "invert" Epsilon using the useful properties of the matrix in R in point-like */
    for( int iR_loc=0; iR_loc < EpsR.get_nblocks(); ++iR_loc ) {
        output::print("iR: ", iR_loc);

        for( int irow = 0; irow < EpsR.get_nrows(); ++irow ) {
            for( int icol = 0; icol < EpsR.get_ncols(); ++icol ) {
                invEpsR(iR_loc, irow, icol) = (std::abs(EpsR(iR_loc, irow, icol)) < 1.e-03) ? 0 : 1./(EpsR(iR_loc, irow, icol));              
            }
        }
    }

    InverseEpsilon_.print_Rdecay("InverseEpsilon", material_.rwann_);
#ifdef EDUS_HDF5
    InverseEpsilon_.go_to_k(true);
    HDF5_tree fout("output.h5", hdf5_access_t::read_write);
    mpi::Communicator::world().barrier();
    InverseEpsilon_.get_Operator_k().write_h5("output.h5",nodename_,"Epsilonk");
    InverseEpsilon_.go_to_R(false);
#endif

    /* get screened interaction W=eps^-1 V in R space. We know how they behave */
    W_.go_to_R(false);
    auto& WR = W_.get_Operator(Space::R);
    auto& VR = V_.get_Operator(Space::R);

    for( int iR_loc=0; iR_loc < EpsR.get_nblocks(); ++iR_loc ) {
        for( int irow = 0; irow < EpsR.get_nrows(); ++irow ) {
            for( int icol = 0; icol < EpsR.get_ncols(); ++icol ) {
                WR(iR_loc, irow, icol) = invEpsR(iR_loc, irow, icol)*VR(iR_loc, irow, icol);
            }
        }
    }

    W_.print_Rdecay("W", material_.rwann_);
#ifdef EDUS_HDF5
    W_.go_to_k(true);
    mpi::Communicator::world().barrier();
    W_.get_Operator_k().write_h5("output.h5",nodename_,"Wk");
    W_.go_to_R(false);
#endif

}
