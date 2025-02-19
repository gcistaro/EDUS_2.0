#include "Screened_potential.hpp"

void ScreenedPotential::initialize()
{
    /* Initialize all the dimensions */
    W_      .initialize_fft(DensityMatrix_);
    V_      .initialize_fft(DensityMatrix_);
    X_      .initialize_fft(DensityMatrix_);
    Epsilon_.initialize_fft(DensityMatrix_);
}

void ScreenedPotential::ResponseFunction()
{
    X_.go_to_k(false);

    /* Lindhard formula in point like approximation to get X_{ij}(q) */
    auto num_bands = Epsilon_.get_Operator(k).get_nrows();
    auto num_filledbands =  ctx_->cfg().filledbands();
    auto num_emptybands  =  num_bands - num_filledbands;
    auto kmesh     = Epsilon_.get_Operator(k).get_MeshGrid();

    static Matrix<std::complex<double>> UdkpUk(num_emptybands, num_filledbands); //contains U^dagger(k')U(k)
    static Matrix<std::complex<double>> UdkUkp(num_filledbands, num_emptybands); //contains U^dagger(k)U(k')
    static Matrix<std::complex<double>> temp(num_emptybands, num_emptybands); //contains the multiplication 
                                                                                        
    for( int iq = 0; iq < kmesh->get_TotalSize(); ++iq ) {
        output::print("Calculation for iq = ", iq);
        auto q = (*kmesh)[iq];
        for ( int ik = 0; ik < kmesh->get_TotalSize(); ++ik ) {
            auto k = (*kmesh)[ik];
            auto kq_reduced = kmesh->reduce(k+q);
            auto ikp = kmesh->find( kq_reduced );//index of k+q

            auto& Uk   = Operator<std::complex<double>>::EigenVectors[ik];
            auto& Udkp = Operator<std::complex<double>>::EigenVectors_dagger[ikp]; // WARNING we need e^{iGr}, but it cancels
            auto& Ukp  = Operator<std::complex<double>>::EigenVectors[ikp]; // WARNING we need e^{-iGr}, but it cancels
            auto& Udk  = Operator<std::complex<double>>::EigenVectors_dagger[ik];

            /* U^dagger(k').U(k) */
            Matrix_gemm(UdkpUk, 1., Udkp, Uk, 0.);

            /* U^\dagger(k).U(k') */
            Matrix_gemm(UdkUkp, 1., Udk, Ukp, 0.); 
            /* UdkUkp_vc /= ( Ec(kp) - Ev(k) )*/
            for ( int ival = 0; ival < num_filledbands; ++ival ) {
                for ( int icond = 0; icond < num_emptybands; ++icond ) {
                    UdkUkp(ival, icond) /= ( Band_energies_[ikp](num_filledbands+icond) - Band_energies_[ik](ival) );
                }
            }

            Matrix_gemm(temp, 1., UdkpUk, UdkUkp, 0.);
        
            //Trace over empty bands
            std::complex<double> Trace = 0.;
            for ( int icond = 0; icond < num_emptybands; ++icond ) {
                Trace += temp(icond, icond);
            }

            for ( int iwann=0; iwann<num_bands; ++iwann ) {
                X_.get_Operator_k()( iq, iwann, iwann ) = 2.*Trace/double(kmesh->get_TotalSize()*kmesh->get_TotalSize());
            }
        }
    }

    X_.go_to_R();
    std::stringstream rank;

#ifdef EDUS_MPI
    rank << "X" << mpi::Communicator::world().rank() << ".txt";
#else
    rank << "X0.txt";
#endif
    std::ofstream os;
    os.open(rank.str());
    auto Rgamma_centered = get_GammaCentered_grid(*X_.get_Operator_R().get_MeshGrid());
    for(int iR_loc=0; iR_loc< DensityMatrix_.get_Operator_R().get_nblocks(); ++iR_loc){
        //os << DensityMatrix.get_Operator_R()[i] << std::endl;
        auto iR_glob = DensityMatrix_.mpindex.loc1D_to_glob1D(iR_loc);
        os << Rgamma_centered[iR_glob].norm();
        os << " " << std::abs(max(DensityMatrix_.get_Operator_R()[iR_loc])) << std::endl;
    }
    os.close();

}

void ScreenedPotential::BareCoulomb()
{
    V_.go_to_R();
    V_.get_Operator(Space::k).fill(0.);

    auto index_origin = material_.r[0].get_Operator_R().get_MeshGrid()->find(Coordinate(0,0,0));
    auto& x0 = (material_.r[0].get_Operator_R())[index_origin];
    auto& y0 = (material_.r[1].get_Operator_R())[index_origin];
    auto& z0 = (material_.r[2].get_Operator_R())[index_origin];

    for( int iR_loc=0; iR_loc < V_.get_Operator_R().get_nblocks(); ++iR_loc ) {
        for( int irow = 0; irow < V_.get_Operator_R().get_nrows(); ++irow ) {
            auto iR_glob = iR_loc;//TODO change for mpi
            auto& Rcart = (*V_.get_Operator_R().get_MeshGrid())[iR_glob].get("Cartesian");

            auto ratom = Coordinate(  x0(irow, irow).real() - Rcart[0],
                                      y0(irow, irow).real() - Rcart[1],
                                      z0(irow, irow).real() - Rcart[2]);
            V_.get_Operator_R()(iR_loc, irow, irow) = 1./ratom.norm();
        }
    }
    V_.go_to_k();
}

void ScreenedPotential::Epsilon()
{
    /* Implementation of Epsilon(r,r') = \delta(r-r') -\int d^3 r'' V(r,r'')X(r'',r') */

    /* -VX <- Epsilon in q space */
    V_.go_to_k();
    multiply(Epsilon_.get_Operator_k(), -1.+0.*im, V_.get_Operator_k(), X_.get_Operator_k());
    
    /* add delta(r-r'). Better to do it in R space there it is delta_{R0} delta_{ij} */
    Epsilon_.go_to_R();
    auto index_origin = Epsilon_.get_Operator_R().get_MeshGrid()->find(Coordinate(0.,0.,0.));
    
    for( int irow = 0; irow < V_.get_Operator_R().get_nrows(); ++irow ) {
        Epsilon_.get_Operator_R()(index_origin, irow, irow) += 1.;
    }
}

void ScreenedPotential::Calculate()
{
    ResponseFunction();
    BareCoulomb();
    Epsilon();

    /* get inverse of epsilon */
    Epsilon_.go_to_k();
    for( int ik_loc=0; ik_loc < Epsilon_.get_Operator_R().get_nblocks(); ++ik_loc ) {
        InverseEpsilon_.get_Operator_k()[ik_loc] = Epsilon_.get_Operator_k()[ik_loc].inverse(); 
    }

    /* get screened interaction W=eps^-1 V in q space */
    multiply(W_.get_Operator_k(), 1.+0.*im, InverseEpsilon_.get_Operator_k(), V_.get_Operator_k());
    W_.go_to_R();
}
