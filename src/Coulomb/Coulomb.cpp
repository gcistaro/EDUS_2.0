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

                                                // BARE COULOMB INTERACTION //
 
    auto BareCoulomb_TB = mdarray<std::complex<double>,3> ( { int( size_MG_global ), nbnd, nbnd } );                                                
 
    // locate and open bare coulomb file
    std::filesystem::path barecoulomb_file_path = std::filesystem::current_path() / "barecoulomb.txt";
    auto barecoulomb_file = ReadFile(barecoulomb_file_path.string());

    // read from the kcw file the R vectors where the Hamiltonian is computed
    std::vector<Coordinate> RkcwGrid;
    std::array<double,3> Rkcw;
    for (int iline = 0; iline < barecoulomb_file.size(); iline++)
    {
        if (barecoulomb_file[iline].size() == 3)
        {
            for (int i = 0; i < 3; i++){
              Rkcw[i] = stoi(barecoulomb_file[iline][i]);}
            Coordinate R(Rkcw[0], Rkcw[1], Rkcw[2], LatticeVectors(Space::R));
            RkcwGrid.push_back(R);
        }
    }

    // find in the systems grid the R vectors from the kcw file read above
    RCoulomb.initialize(Space::R, RkcwGrid, 0.0);
    auto Rgrid_shifted = get_GammaCentered_grid(*Rgrid);
    MeshGrid::Calculate_ConvolutionIndex(RCoulomb, Rgrid_shifted, *Operator<std::complex<double>>::MeshGrid_Null);
    auto& ci = MeshGrid::ConvolutionIndex[{RCoulomb.get_id(), Rgrid_shifted.get_id(), Operator<std::complex<double>>::MeshGrid_Null->get_id()}];

    // build the bare coulomb interaction matrix elements in the imported R vectors
    BareCoulomb_TB.fill(0.0);
    for (int iRCoulomb=0; iRCoulomb<RCoulomb.get_TotalSize(); iRCoulomb++)
    {
        for (int irow=0; irow<nbnd; irow++)
        {
            for (int icol=0; icol<nbnd; icol++)
            {
                int iline = nbnd*2*irow + 2*icol + (std::pow(nbnd,2)*2+1)*iRCoulomb + 1; // the 2 is due to spin

                BareCoulomb_TB(ci(iRCoulomb,0), irow, icol) = Convert(std::atof(barecoulomb_file[iline][3].c_str()) + 
                std::atof(barecoulomb_file[iline+1][3].c_str()) - 0.14891086, Hartree, AuEnergy); // the constant is due to the way Nicola computes something ¯\_(ツ)_/¯
                //std::cout << "BareCoulomb_TB(ci(" << iRCoulomb << ",0), " << irow << ", " << icol << ") = " << BareCoulomb_TB(ci(iRCoulomb,0), irow, icol) << std::endl;
            }
        }
    }

    Operator<std::complex<double>> Bare;
    Bare.initialize_fft(*Rgrid__, nbnd);
    auto& bareR = Bare.get_Operator(Space::R);
    for (int iRCoulomb=0; iRCoulomb<RCoulomb.get_TotalSize(); iRCoulomb++)// maybe Rgrid__
    {
        for (int irow=0; irow<nbnd; irow++)
        {
            for (int icol=0; icol<nbnd; icol++)
            {
                bareR(iRCoulomb, irow, icol ) = BareCoulomb_TB(iRCoulomb, irow, icol);
            }
        }
    }

    std::vector<Coordinate> rwann(nbnd);
    for (auto& rwann_iwann : rwann) {
        rwann_iwann.initialize(0.,0.,0.);
    }
    Bare.print_Rdecay("EDUSbarecoulomb",rwann);

                                            // SCREENED COULOMB INTERACTION //
    
    auto ScreenCoulomb_TB = mdarray<std::complex<double>,3> ( { int( size_MG_global ), nbnd, nbnd } );

    // import screened coulomb interaction
    std::filesystem::path screencoulomb_file_path = std::filesystem::current_path() / "screencoulomb.txt";
    auto screencoulomb_file = ReadFile(screencoulomb_file_path.string());

    // build the screened coulomb interaction matrix elements in the imported R vectors
    ScreenCoulomb_TB.fill(0.0);
    for (int iRCoulomb=0; iRCoulomb<RCoulomb.get_TotalSize(); iRCoulomb++)
    {
        for (int irow=0; irow<nbnd; irow++)
        {
            for (int icol=0; icol<nbnd; icol++)
            {
                int iline = nbnd*2*irow + 2*icol + (std::pow(nbnd,2)*2+1)*iRCoulomb + 1; // the two is due to spin

                ScreenCoulomb_TB(ci(iRCoulomb,0), irow, icol) = Convert(std::atof(screencoulomb_file[iline][3].c_str()) + 
                std::atof(screencoulomb_file[iline+1][3].c_str()), Hartree, AuEnergy);
                //std::cout << "ScreenCoulomb_TB(ci(" << iRCoulomb << ",0), " << irow << ", " << icol << ") = " << ScreenCoulomb_TB(ci(iRCoulomb,0), irow, icol) << std::endl;
            }
        }
    }

    Operator<std::complex<double>> Screen;
    Screen.initialize_fft(*Rgrid__, nbnd);
    auto& screenR = Screen.get_Operator(Space::R);
    for (int iRCoulomb=0; iRCoulomb<RCoulomb.get_TotalSize(); iRCoulomb++)// maybe Rgrid__
    {
        for (int irow=0; irow<nbnd; irow++)
        {
            for (int icol=0; icol<nbnd; icol++)
            {
                screenR(iRCoulomb, irow, icol ) = ScreenCoulomb_TB(iRCoulomb, irow, icol);
            }
        }
    }

    Screen.print_Rdecay("EDUSscreencoulomb",rwann);



                                                    // COMBINING THE TWO TERMS //

    /* Get local part and add the minus sign */
    #pragma omp parallel for
    for( int iR_local = 0; iR_local < size_MG_local; ++iR_local )
    {
        for( int irow = 0; irow < nbnd; ++irow )
        {
            for( int icol = 0; icol < nbnd; ++icol )
            {
                auto iR_global = int( Rgrid__->mpindex.loc1D_to_glob1D(iR_local) );
                HF( iR_local, irow, icol ) = BareCoulomb_TB( iR_global, irow, icol ) - ScreenCoulomb_TB( iR_global, irow, icol ); 
            }
        }
    }
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
