

#include "Coulomb.hpp"
#include "RytovaKeldysh/RytovaKeldysh.hpp"

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

void Coulomb::initialize(const int& nbnd, const std::shared_ptr<MeshGrid>& Rgrid__, const std::array<Operator<std::complex<double>>,3>& r)
{
    auto size_MG =  Rgrid__->get_TotalSize();
    //W = mdarray<std::complex<double>, 6>( { size_MG, nbnd, nbnd, size_MG, nbnd, nbnd } );
    HF = mdarray<std::complex<double>,3> ( { int( size_MG ), nbnd, nbnd } );
    /* initializing W in point like approximation */ 
    RytovaKeldysh RytKel;
    RytKel.initialize(r, 2, Rgrid__);
    //Fock part
    #pragma omp parallel for
    for( int iR = 0; iR < size_MG; ++iR ){
        for( int irow = 0; irow < nbnd; ++irow ){
            for( int icol = 0; icol < nbnd; ++icol ){
                HF( iR, irow, icol ) = - RytKel.TB( iR, irow, icol ); 
            }
        }
    }

    int index_origin = Rgrid__->find(Coordinate(0,0,0));
    std::cout << index_origin;
    #pragma omp parallel for
    for( int irow = 0; irow < nbnd; ++irow ){
        HF( index_origin, irow, irow ) += 2.*RytKel.TB( index_origin, irow, irow );
    }
    std::cout << "Done HF!\n";
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
    if ( !DoCoulomb ) {
        return;
    }

    auto& HR = H__.get_Operator(R);
    auto& DMR = DM__.get_Operator(R);
    auto& DM0R = DM0.get_Operator(R);

    assert( HR.get_MeshGrid()->get_id() == DMR.get_MeshGrid()->get_id() );

    if( EraseH__ ) {
        HR.fill(0.);
    }

    //PointLike approximation
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



