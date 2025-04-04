#include "Coulomb.hpp"
#include "ModelCoulomb/ModelCoulomb.hpp"
#include <filesystem> 

/// @brief Triggers the initialization of the object
/// @param nbnd__ Number of bands/wannier functions used in the simulation
/// @param Rgrid__ Grid in R space used in the simulation
/// @param r__ Position operator read from Wannier90, in a.u.
Coulomb::Coulomb(const int& nbnd__, const std::shared_ptr<MeshGrid>& Rgrid__, const std::array<Operator<std::complex<double>>,3>& r__)
{
    initialize(nbnd__, Rgrid__, r__);
}
/// @brief Initialize the objects of the class, mainly ModelCoulomb
/// @param nbnd__ Number of bands/wannier functions used in the simulation
/// @param Rgrid__ Grid in R space used in the simulation
/// @param r__ Position operator read from Wannier90, in a.u.
void Coulomb::initialize(const int& nbnd, const std::shared_ptr<MeshGrid>& Rgrid__, const std::array<Operator<std::complex<double>>,3>& r__)
{
    if(!DoCoulomb_) {
        return;
    }
    Rgrid_ = Rgrid__;

    // == auto size_MG_global =  Rgrid__->get_TotalSize();
    // == auto size_MG_local = Rgrid__->get_LocalSize();
    // == HF = mdarray<std::complex<double>,3> ( { int( size_MG_local ), nbnd, nbnd } );
    // == std::filesystem::path cwd = std::filesystem::current_path() / "RytovaKeldysh.txt";
    // == read_rk_py( RytovaKeldysh_TB, cwd.str());

    modelcoulomb_.initialize(r__, 2, Rgrid__);
    int index_origin_global = Rgrid__->find(Coordinate(0,0,0));
    HasOrigin_ = Rgrid__->mpindex.is_local(index_origin_global);
    if( HasOrigin_ ) {
        index_origin_local_ = Rgrid__->mpindex.glob1D_to_loc1D(index_origin_global);  
    }  
}

/// @brief Setter for DM0 (Density Matrix of the ground state at Wannier gauge in R)
/// @param DM0__ Value to use in the setter
void Coulomb::set_DM0( const Operator<std::complex<double>>& DM0__ )
{
    DM0_ = DM0__;
}

/// @brief Setter for Docoulomb variable of the class 
/// @param DoCoulomb__ Value to use in the setter
void Coulomb::set_DoCoulomb(const bool& DoCoulomb__)
{
    DoCoulomb_ = DoCoulomb__;
}

/// @brief Set Epsilon (macroscopic dielectric constant)
/// @param Epsilon__ Value we want to use as epsilon
void Coulomb::set_epsilon(const bool& Epsilon__)
{
    modelcoulomb_.set_epsilon( Epsilon__ );
}

/// @brief Set r0 in RytovaKeldysh model, not used otherwise
/// @param r0__ Value we want to use as r0 (in a.u.)
void Coulomb::set_r0(const bool& r0__)
{
    modelcoulomb_.set_r0( r0__ );
}

/// @brief Getter for DoCoulomb variable 
/// @return DoCoulomb variable
const bool& Coulomb::get_DoCoulomb() const
{
    return DoCoulomb_;
}

/// @brief Getter for DoCoulomb variable 
/// @return DoCoulomb variable
bool& Coulomb::get_DoCoulomb()
{
    return DoCoulomb_;
}

/// @brief This function calculates the effective Hamiltonian from the Coulomb interaction. 
/// The Coulomb interaction has two different terms: 
/// - The Hartree term 
/// @f[
///   H^{(\text{eff})}_{n'n}(\textbf{R}) = \delta_{nn'} \delta_{\bf R, 0} \sum_{m\bf R'} V_{nm}[\textbf{R}']\rho_{mm}(\textbf{0}) 
/// @f]
/// - The Fock term
/// @f[
///   H^{(\text{eff})}_{n'n}(\textbf{R}) = -W_{n'n}[\textbf{R}]\rho_{n'n}(\textbf{R}) 
/// @f]
/// @param H__ Hamiltonian Operator, in which we want to add the effective Hamiltonian of the Coulomb interaction
/// @param DM__ Density matrix at current time, that we need to use to calculate the effective Coulomb interaction
/// @param EraseH__ True if we want to erase the H__ matrix before feeding it with the effective Coulomb interaction
void Coulomb::EffectiveHamiltonian(Operator<std::complex<double>>& H__, const Operator<std::complex<double>>& DM__,
                                  const bool& EraseH__ ) 
{
    /* We calculate the Coulomb interaction in R space */
    auto& HR__ = H__.get_Operator(R);
    auto& DMR__ = DM__.get_Operator(R);
    auto& DM0R_ = DM0_.get_Operator(R);

    if( EraseH__ ) {
        HR__.fill(0.);
    }

    if ( !DoCoulomb_ ) {
        return;
    }

    /* Hartree term */
    if( HasOrigin_ ) { // Only the rank with R=0 contributes to this term 
        #pragma omp parallel for
        for( int iblock = 0; iblock < HR__.get_nblocks(); ++iblock ) {
            for( int irow = 0; irow < HR__.get_nrows(); ++irow ) {
                for( int icol = 0; icol < HR__.get_ncols(); ++icol ) {
                    HR__(index_origin_local_, irow, irow) += 
                            modelcoulomb_.BarePotential_(iblock, irow, icol)*(DMR__(iblock, icol, icol) - DM0R_(iblock, icol, icol)); 
                }
            }
        }
    }

    /* Fock term */
    #pragma omp parallel for
    for( int iblock = 0; iblock < HR__.get_nblocks(); ++iblock ) {
        for( int irow = 0; irow < HR__.get_nrows(); ++irow ) {
            for( int icol = 0; icol < HR__.get_ncols(); ++icol ) {
                HR__( iblock, irow, icol ) -= 
                        modelcoulomb_.ScreenedPotential_( iblock, irow, icol )*( DMR__( iblock, irow, icol ) - DM0R_( iblock, irow, icol ) );
            }
        }
    } 

    //for( int iblock = 0; iblock < H_.get_nblocks(); ++iblock ) {
    //    for( int irow = 0; irow < H_.get_nrows(); ++irow ) {
    //        for( int icol = 0; icol < H_.get_ncols(); ++icol ) {
    //            H_( iblock, irow, icol ) += W( iblock, irow, icol, iblock1, ibnd1, ibnd2 )*DM( iblock1, ibnd1, ibnd2 );
    //        }
    //    }
    //} 
}    


void read_rk_py(mdarray<std::complex<double>,3>& RytovaKeldysh_TB, const std::string& filename)
{
    auto file = ReadFile(filename);
    auto index = 1;
    
    /* Read Rytova Keldysh (screened) potential produced with RytovaKeldysh.py */
    for( int iR = 0; iR < int(RytovaKeldysh_TB.get_Size(0)); ++iR ) {
        for( int irow = 0; irow < int(RytovaKeldysh_TB.get_Size(1)); ++irow ) {
            for( int icol = 0; icol < int(RytovaKeldysh_TB.get_Size(1)); ++icol ) {
                assert(file[index].size() == 1);
                RytovaKeldysh_TB( iR, irow, icol ) = std::atof( file[index][0].c_str() );
                index++;
            }
        }
    }
    assert(index == int(file.size()));
}
