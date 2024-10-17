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

void Coulomb::initialize(const int& nbnd, const std::shared_ptr<MeshGrid>& Rgrid__, const std::array<Operator<std::complex<double>>,3>& r)
{
    if(!DoCoulomb) {
        return;
    }
    Rgrid = Rgrid__;
    auto size_MG_global =  Rgrid__->get_TotalSize();
    auto size_MG_local = Rgrid__->get_LocalSize();
    //W = mdarray<std::complex<double>, 6>( { size_MG, nbnd, nbnd, size_MG, nbnd, nbnd } );
    HF = mdarray<std::complex<double>,3> ( { int( size_MG_local ), nbnd, nbnd } );
    /* initializing W in point like approximation */ 
    auto RytovaKeldysh_TB = mdarray<std::complex<double>,3> ( { int( size_MG_global ), nbnd, nbnd } );
/*
    RytovaKeldysh RytKel;
    RytKel.initialize(r, 2, Rgrid__);

*/
    std::filesystem::path cwd = std::filesystem::current_path() / "RytovaKeldysh.txt";
    auto file = ReadFile(cwd.string());
    auto index = 1;
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

    //Fock part
    #pragma omp parallel for
    for( int iR_local = 0; iR_local < size_MG_local; ++iR_local ){
        for( int irow = 0; irow < nbnd; ++irow ){
            for( int icol = 0; icol < nbnd; ++icol ){
                auto iR_global = int( Rgrid__->mpindex.loc1D_to_glob1D(iR_local) );
                HF( iR_local, irow, icol ) = - RytovaKeldysh_TB( iR_global, irow, icol ); 
            }
        }
    }

    int index_origin_global = Rgrid__->find(Coordinate(0,0,0));
    if( Rgrid__->mpindex.is_local(index_origin_global) ) {
        int index_origin_local = Rgrid__->mpindex.glob1D_to_loc1D(index_origin_global);

        for( int iR = 0; iR < int(RytovaKeldysh_TB.get_Size(0)); ++iR ) {
            for( int irow = 0; irow < nbnd; ++irow ){
                HF( index_origin_local, irow, irow ) += 2.*RytovaKeldysh_TB( iR, irow, irow );
            }
        }
    }
        std::ofstream HFF("HF.txt");
    for( int iR = 0; iR < size_MG_local; ++iR ){
        for( int irow = 0; irow < nbnd; ++irow ){
            for( int icol = 0; icol < nbnd; ++icol ){
                HFF << HF( iR, irow, icol ).real() << " " << HF( iR, irow, icol ).imag() <<  std::endl;
            }}}
            HFF.close();

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



