

#include "Coulomb.hpp"

Coulomb::Coulomb(const int& nbnd, const MeshGrid& Rgrid__)
{
    initialize(nbnd, Rgrid__);
}


void Coulomb::initialize(const int& nbnd, const MeshGrid& Rgrid__)
{
    auto size_MG = Rgrid__.get_TotalSize();
    W = mdarray<std::complex<double>, 6>( { size_MG, nbnd, nbnd, size_MG, nbnd, nbnd } );

    /* initializing W */ 
}

void Coulomb::EffectiveHamiltonian(BlockMatrix<std::complex<double>>& H_, const BlockMatrix<std::complex<double>>& DM,
                                  const bool& EraseH_ ) 
{
    assert( H_.get_MeshGrid()->get_id() == DM.get_MeshGrid()->get_id() == this->get_MeshGrid() );

    if( EraseH_ ) {
        H_.fill(0.);
    }

    for( int iblock = 0; iblock < H_.get_nblocks(); ++iblock ) {
        for( int irow = 0; irow < H_.get_nrows(); ++irow ) {
            for( int icol = 0; icol < H_.get_ncols(); ++icol ) {
                H_( iblock, irow, icol ) += W( iblock, irow, icol, iblock1, ibnd1, ibnd2 )*DM( iblock1, ibnd1, ibnd2 );
            }
        }
    } 
}    



