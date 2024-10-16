

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
    if(!DoCoulomb) {
        return;
    }
    auto size_MG =  Rgrid__->get_TotalSize();
    //W = mdarray<std::complex<double>, 6>( { size_MG, nbnd, nbnd, size_MG, nbnd, nbnd } );
    HF = mdarray<std::complex<double>,3> ( { int( size_MG ), nbnd, nbnd } );
    /* initializing W in point like approximation */ 
    auto RytovaKeldysh_TB = mdarray<std::complex<double>,3> ( { int( size_MG ), nbnd, nbnd } );
/*
    RytovaKeldysh RytKel;
    RytKel.initialize(r, 2, Rgrid__);

*/
    auto det = Coordinate::get_Basis(LatticeVectors(R)).get_M().determinant()/Coordinate::get_Basis(LatticeVectors(R)).get_M()(2,2);
    std::ifstream("RytovaKeldysh.txt");
    auto file = ReadFile("RytovaKeldysh.txt");
std::cout << "Coordinate::get_Basis(LatticeVectors(R)).get_M().determinant():" << det << " " << std::sqrt(det) << std::endl;//.determinant() << std::endl;
    auto index = 1;
    for( int iR = 0; iR < int(HF.get_Size(0)); ++iR ) {
        for( int irow = 0; irow < int(HF.get_Size(1)); ++irow ) {
            for( int icol = 0; icol < int(HF.get_Size(1)); ++icol ) {
                std::cout << "index: " << index << std::endl;
                assert(file[index].size() == 1);
                RytovaKeldysh_TB( iR, irow, icol ) = std::atof( file[index][0].c_str() );
                //RytovaKeldysh_TB( iR, irow, icol ) *= 100.*100.*100.*100.;
                index++;
            }
        }
    }
    assert(index == int(file.size()));

    //Fock part
    #pragma omp parallel for
    for( int iR = 0; iR < size_MG; ++iR ){
        for( int irow = 0; irow < nbnd; ++irow ){
            for( int icol = 0; icol < nbnd; ++icol ){
                HF( iR, irow, icol ) = - RytovaKeldysh_TB( iR, irow, icol ); 
            }
        }
    }
    int index_origin = Rgrid__->find(Coordinate(0,0,0));
    std::cout << index_origin;
    #pragma omp parallel for
    for( int iR = 0; iR < int(HF.get_Size(0)); ++iR ) {
        for( int irow = 0; irow < nbnd; ++irow ){
            HF( index_origin, irow, irow ) += 2.*RytovaKeldysh_TB( iR, irow, irow );
        }
    }
    std::cout << HF;
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
    auto& HR = H__.get_Operator(R);
    auto& DMR = DM__.get_Operator(R);
    auto& DM0R = DM0.get_Operator(R);

    if( EraseH__ ) {
        HR.fill(0.);
    }

    if ( !DoCoulomb ) {
        return;
    }
BlockMatrix<std::complex<double>> Delta_DM;
Delta_DM = DMR;
Delta_DM.fill(0);
    #pragma omp parallel for
    for( int iblock = 0; iblock < HR.get_nblocks(); ++iblock ) {
        for( int irow = 0; irow < HR.get_nrows(); ++irow ) {
            for( int icol = 0; icol < HR.get_ncols(); ++icol ) {
                HR( iblock, irow, icol ) += HF( iblock, irow, icol )*( DMR( iblock, irow, icol ) - DM0R( iblock, irow, icol ) );
            }
        }
    } 
    //std::cout << "max(Delta_DM): " << *max(Delta_DM)<<std::endl;
    /*
    auto& HR = H__.get_Operator(R);
    auto& DMR = DM__.get_Operator(R);
    auto& DM0R = DM0.get_Operator(R);

    if( EraseH__ ) {
        HR.fill(0.);
    }

    if ( !DoCoulomb ) {
        return;
    }

    assert( HR.get_MeshGrid()->get_id() == DMR.get_MeshGrid()->get_id() );

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
*/
}    



