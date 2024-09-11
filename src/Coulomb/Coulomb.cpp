

#include "Coulomb.hpp"

Coulomb::Coulomb(const int& nbnd, const MeshGrid& Rgrid__)
{
    initialize(nbnd, Rgrid__);
}


void Coulomb::initialize(const int& nbnd, const std::shared_ptr<MeshGrid>& Rgrid__)
{
    Rgrid = Rgrid__;
    auto size_MG = Rgrid__.get_TotalSize();
    HF = mdarray<std::complex<double>, 6>( { size_MG, nbnd, nbnd, size_MG, nbnd, nbnd } );

    /* initializing W and W0*/ 
    auto Wtot = mdarray<std::complex<double>,7>({size_MG, size_MG, size_MG, nbnd, nbnd, nbnd, nbnd});
    
    RytovaKeldysh( 2, Rgrid);
    
    Read_All_W(Wtot);//get full Coulomb potential: 4 band indices, 3 R indices
    HartreeFock(Wtot, HF);//get Hartree and Fock potential from Wtot
    

}

void Coulomb::EffectiveHamiltonian(BlockMatrix<std::complex<double>>& H_, const BlockMatrix<std::complex<double>>& DM,
                                  const bool& EraseH_ ) 
{
    assert( H_.get_MeshGrid()->get_id() == DM.get_MeshGrid()->get_id() == this->get_MeshGrid() );

    if( EraseH_ ) {
        H_.fill(0.);
    }

    #pragma omp parallel for schedule(static)
    for( int iblock = 0; iblock < H_.get_nblocks(); ++iblock ) {
        for( int irow = 0; irow < H_.get_nrows(); ++irow ) {
            for( int icol = 0; icol < H_.get_ncols(); ++icol ) {
                for( int iblock1 = 0; iblock1 < H_.get_nblocks(); ++iblock1 ) {
                    for( int irow1 = 0; irow1 < H_.get_nrows(); ++irow1 ) {
                        for( int icol1 = 0; icol1 < H_.get_ncols(); ++icol1 ) {
                            H_( iblock, irow, icol ) +=
                            //this is multiplication of matrix and vector (group indices 3)
                                HF( iblock, irow, icol, iblock1, irow1, icol1 )*DM( iblock1, irow1, icol1 );
                        }
                    }
                }
            }
        }
    } 
}    


//this function is called only once at the beginning
void HartreeFock(const mdarray<std::complex<double>,7>& Wtot, mdarray<std::complex<double>,6> W)
{
    HF.fill(0.);
    #pragma omp parallel for schedule(static)
    for( int iR = 0; iR < H_.get_nblocks(); ++iR ) {
        for( int in = 0; in < H_.get_nrows(); ++in ) {
            for( int in_ = 0; in_ < H_.get_ncols(); ++in_ ) {
                for( int iT = 0; iT < H_.get_nblocks(); ++iT ) {
                    for( int im = 0; im < H_.get_nrows(); ++im ) {
                        for( int im_ = 0; im_ < H_.get_ncols(); ++im_ ) {
                            for( int iS = 0; iS < H_.get_nblocks(); ++iS ) {
                              HF(iR,in,in_,iT,im,im_) +=
                                 2*Wtot(in,ci(iR,iS),im,iT,in_,iS,im_) //Hartree
                                -Wtot(in,iR,im,ci(iT,iS),im_,iS,in_);  //Fock
                            }
                        }
                    }
                }
            }
        }
    }
}


     auto Wtot = mdarray<std::complex<double>,7>({size_MG, size_MG, size_MG, nbnd, nbnd, nbnd, nbnd});
void Read_All_W(mdarray<std::complex<double>,7>& Wtot)
{
    Wtot.fill(0.);

    #pragma omp parallel for schedule(static)
    for( int iR = 0; iR < H_.get_nblocks(); ++iR ) {
        for( int in = 0; in < H_.get_nrows(); ++in ) {
            for( int im = 0; im < H_.get_nrows(); ++im ) {
                Wtot(iR,iR,0,in,in,im,im) = RytovaKeldysh.TB(iR,0,in,im); 
            }
        }
    }
}
