#include "Coulomb.hpp"
#include "RytovaKeldysh/RytovaKeldysh.hpp"
#include <filesystem> 

Coulomb::Coulomb(const int& nbnd, const std::shared_ptr<MeshGrid>& Rgrid__)
{
    initialize(nbnd, Rgrid__);
}

Coulomb::Coulomb(const int& nbnd, const std::shared_ptr<MeshGrid>& Rgrid__, const std::array<Operator<std::complex<double>>,3>& r)
{
    initialize(nbnd, Rgrid__, r);
}

void Coulomb::initialize(const int& nbnd, const std::shared_ptr<MeshGrid>& Rgrid__)
{}

void Coulomb::initialize(const int& nbnd, const std::shared_ptr<MeshGrid>& Rgrid__, const std::array<Operator<std::complex<double>>,3>& r__)
{
    if(!DoCoulomb) {
        return;
    }
    Rgrid = Rgrid__;
    auto size_MG_global =  Rgrid__->get_TotalSize();
    auto size_MG_local = Rgrid__->get_LocalSize();
    HF = mdarray<std::complex<double>,3> ( { int( size_MG_local ), nbnd, nbnd } );

    //old version -> rytovekeldysh read from python output
    //auto RytovaKeldysh_TB = mdarray<std::complex<double>,3> ( { int( size_MG_global ), nbnd, nbnd } );
    //std::filesystem::path cwd = std::filesystem::current_path() / "RytovaKeldysh.txt";
    //read_rk_py( RytovaKeldysh_TB, cwd.str());

    rytovakeldysh.initialize(r__, 3, Rgrid__);
    auto& RytovaKeldysh_TB = rytovakeldysh.TB;

    /* Get local part and add the minus sign */
    #pragma omp parallel for
    for( int iR_local = 0; iR_local < size_MG_local; ++iR_local ){
        for( int irow = 0; irow < nbnd; ++irow ){
            for( int icol = 0; icol < nbnd; ++icol ){
                auto iR_global = int( Rgrid__->mpindex.loc1D_to_glob1D(iR_local) );
                HF( iR_local, irow, icol ) = -RytovaKeldysh_TB( iR_global, irow, icol );
            }
        }
    }

    Operator<std::complex<double>> Screen;
     Screen.initialize_fft(*Rgrid__, nbnd);
     auto& screenR = Screen.get_Operator(Space::R);
     for (int iRCoulomb=0; iRCoulomb<(*Rgrid__).get_TotalSize(); iRCoulomb++)// maybe Rgrid__
     {
         for (int irow=0; irow<nbnd; irow++)
         {
             for (int icol=0; icol<nbnd; icol++)
             {
                 screenR(iRCoulomb, irow, icol ) = RytovaKeldysh_TB(iRCoulomb, irow, icol);
             }
         }
     }
 
     std::vector<Coordinate> rwann(nbnd);
     for (auto& rwann_iwann : rwann) {
         rwann_iwann.initialize(0.,0.,0.);
     }
     Screen.print_Rdecay("EDUSmodelscreencoulomb",rwann);

    /* Get unscreened (only on-site for point-like approximation) potential */
    auto index_origin_dft = r__[0].get_Operator(R).get_MeshGrid()->find(Coordinate(0,0,0));
    int index_origin_global = Rgrid__->find(Coordinate(0,0,0));
    std::array<double, 3> bare_ratom;
    if( Rgrid__->mpindex.is_local(index_origin_global) ) {
        int index_origin_local = Rgrid__->mpindex.glob1D_to_loc1D(index_origin_global);

        for( int iR = 0; iR < size_MG_global; ++iR ) {
            auto& R = (*Rgrid__)[iR];
            for( int irow = 0; irow < nbnd; ++irow ){
                for( auto& ix : {0,1,2}) {
                    bare_ratom[ix] = real(r__[ix].get_Operator(Space::R)(index_origin_dft,irow,irow));
                }

                auto ratom = Coordinate(bare_ratom[0], bare_ratom[1], bare_ratom[2]);
                auto norm = ratom.norm();

                HF( index_origin_local, irow, irow ) += ( ratom.norm() < 1.e-06 
                                                                        ? 0. 
                                                                        : 1./((ratom + R).norm()) *2.); //2 for spin degeneracy
            }
        }
    }

    Operator<std::complex<double>> Bare;
     Bare.initialize_fft(*Rgrid__, nbnd);
     auto& bareR = Bare.get_Operator(Space::R);
     for (int iRCoulomb=0; iRCoulomb<(*Rgrid__).get_TotalSize(); iRCoulomb++)// maybe Rgrid__
     {
         for (int irow=0; irow<nbnd; irow++)
         {
             for (int icol=0; icol<nbnd; icol++)
             {
                 bareR(iRCoulomb, irow, icol ) = HF(iRCoulomb, irow, icol)+RytovaKeldysh_TB(iRCoulomb, irow, icol );
             }
         }
     }
 
     Bare.print_Rdecay("EDUSpointlikebarecoulomb",rwann);
}



void Coulomb::set_DM0( const Operator<std::complex<double>>& DM0__ )
{
    DM0 = DM0__;
}

void Coulomb::set_DoCoulomb(const bool& DoCoulomb__)
{
    DoCoulomb = DoCoulomb__;
}

const bool& Coulomb::get_DoCoulomb() const
{
    return DoCoulomb;
}

bool& Coulomb::get_DoCoulomb()
{
    return DoCoulomb;
}


void Coulomb::EffectiveHamiltonian(Operator<std::complex<double>>& H__, const Operator<std::complex<double>>& DM__,
                                  const bool& EraseH__ ) 
{
    auto& HR = H__.get_Operator(R);
    auto& DMR = DM__.get_Operator(R);
    auto& DM0R = DM0.get_Operator(R);

    if( EraseH__ ) {
        HR.fill(0.);
    }

    if ( !DoCoulomb ) {
        return;
    }

    /* Point-like approximation */
    #pragma omp parallel for
    for( int iblock = 0; iblock < HR.get_nblocks(); ++iblock ) {
        for( int irow = 0; irow < HR.get_nrows(); ++irow ) {
            for( int icol = 0; icol < HR.get_ncols(); ++icol ) {
                HR( iblock, irow, icol ) += HF( iblock, irow, icol )*( DMR( iblock, irow, icol ) - DM0R( iblock, irow, icol ) );
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
