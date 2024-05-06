#include "Simulation/Simulation.hpp"


void Simulation::SettingUp_EigenSystem()
{
    PROFILE("Simulation::SetEigensystem");
    //std::cout << "Hamiltonian: " << material.H.get_Operator_R();
    //diagonalizing Hk on the k vectors of the grid of the Density matrix.
    auto& MasterkGrid = DensityMatrix.get_Operator_k().get_MeshGrid()->get_mesh();
    //for(auto& k : MasterkGrid){
    //    std::cout << k.get("LatticeVectors")[0] << " " << k.get("LatticeVectors")[1] << " " << k.get("LatticeVectors")[2] << std::endl;
    //}
    material.H.dft(MasterkGrid, +1);
    auto& Uk = Operator<std::complex<double>>::EigenVectors;
    auto& UkDagger = Operator<std::complex<double>>::EigenVectors_dagger;

    material.H.get_Operator_k().diagonalize(Band_energies, Uk);

    UkDagger.initialize(Uk.get_nblocks(), Uk.get_ncols(), Uk.get_nrows());
    for(int ik=0; ik<UkDagger.get_nblocks(); ++ik){
        for(int ir=0; ir<UkDagger.get_nrows(); ++ir){
            for(int ic=0; ic<UkDagger.get_ncols(); ++ic){
                UkDagger[ik](ir, ic) = conj(Uk[ik](ic,ir));
            }
        }
    }
    //std::cout << "testing U and Udagger...\n";
    //BlockMatrix<std::complex<double>, k> id(Uk.get_nblocks(), Uk.get_nrows(), Uk.get_ncols());
    //multiply(id, std::complex<double>(1.), Uk, UkDagger);
    //for(int ik=0; ik<UkDagger.get_nblocks(); ++ik){
    //    std::cout << ik << " \n"<< id[ik]; 
    //}
};



void Simulation::print_grids()
{
    auto& MasterRgrid = DensityMatrix.get_Operator_R().get_MeshGrid();
    auto& kgrid = DensityMatrix.get_Operator_k().get_MeshGrid();
    auto Rfft = std::make_shared<MeshGrid<R>>(fftPair<k, R>(*kgrid));//Rspace as fourier space of kgrid
    std::ofstream os;
    os.open("MasterR.txt");
    for (auto &R : MasterRgrid->get_mesh())
    {
        os << R.get("Cartesian");
    }
    os.close();
    os.open("k.txt");
    for (auto &k : kgrid->get_mesh())
    {
        os << k.get("Cartesian");
    }
    os.close();
    auto MeshRfft = fftPair<k, R>(*kgrid);
    os.open("R.txt");
    for (auto &R : Rfft->get_mesh())
    {
        os << R.get("Cartesian");
    }
    os.close();
}
        

void Simulation::Calculate_TDHamiltonian(const double& time)
{
    PROFILE("SourceTerm::Calculate_TDHamiltonian");
    auto& HR = H.get_Operator_R();
    auto& H0R = material.H.get_Operator_R();
    auto& xR = material.r[0].get_Operator_R();
    auto& yR = material.r[1].get_Operator_R();
    auto& zR = material.r[2].get_Operator_R();
    auto las = laser(time).get("Cartesian");
    HR.fill(0.);

    auto& ci = MeshGrid<R>::ConvolutionIndex1[{H0R.get_MeshGrid()->get_id(), 
                                              HR.get_MeshGrid()->get_id(), 
                                              Operator<std::complex<double>>::MeshGrid_Null->get_id()}];
    if(ci.get_Size(0) == 0 ){
        MeshGrid<R>::Calculate_ConvolutionIndex1( *(H0R.get_MeshGrid()), *(HR.get_MeshGrid()), *(Operator<std::complex<double>>::MeshGrid_Null));
    }

    for(int iblock=0; iblock<H0R.get_nblocks(); ++iblock){
        for(int irow=0; irow<H0R.get_nrows(); ++irow){
            for(int icol=0; icol<H0R.get_ncols(); ++icol){
                HR(ci(iblock, 0), irow, icol) = H0R(iblock, irow, icol)
                                       + las[0]*xR(iblock, irow, icol)
                                       + las[1]*yR(iblock, irow, icol)
                                       + las[2]*zR(iblock, irow, icol);
            }
        }
    }
}

void Simulation::Propagate()
{
    std::ofstream Pop, Las, os_Pos;
    //Pop.open("Population.txt");
    //Las.open("Laser.txt");
    //os_Pos.open("Position.txt");

        PROFILE("RK_Propagate");
        //DensityMatrix.go_to_k();
        //DensityMatrix.go_to_bloch();
        //
        //std::vector<std::complex<double>> Population(DensityMatrix.get_Operator_k().get_nrows(), 0.);
        //for(int ik=0; ik<DensityMatrix.get_Operator_k().get_nblocks(); ik++){
        //    for(int ibnd=0; ibnd<Population.size(); ibnd++){
        //        Population[ibnd] += DensityMatrix.get_Operator_k()[ik](ibnd,ibnd);
        //    }
        //}
        //Population[0] = 1.-Population[0]/double(DensityMatrix.get_Operator_k().get_MeshGrid()->get_TotalSize());
        //Population[1] /= double(DensityMatrix.get_Operator_k().get_MeshGrid()->get_TotalSize());
        
        //Las << laser(i*RK_object.get_ResolutionTime()).get("Cartesian");

        //for(int ibnd=0; ibnd<Population.size(); ibnd++){
        //    Pop << std::setw(20) << std::setprecision(10) << Population[ibnd].real();
        //    Pop << ' ';
        //}
        //Pop << std::endl;
//
        //std::array<std::complex<double>, 3> v;
        //for(auto ix : {0, 1, 2}){
        //    v[ix] = Trace(Velocity[ix].get_Operator_R(), DensityMatrix.get_Operator_R());
        //}
        //os_Pos << std::setw(20) << std::setprecision(8) << v[0].real();
        //os_Pos << std::setw(20) << std::setprecision(8) << v[0].imag();
        //os_Pos << std::setw(20) << std::setprecision(8) << v[1].real();
        //os_Pos << std::setw(20) << std::setprecision(8) << v[1].imag();
        //os_Pos << std::setw(20) << std::setprecision(8) << v[2].real();
        //os_Pos << std::setw(20) << std::setprecision(8) << v[2].imag();
        //os_Pos << std::endl;

        //DensityMatrix.go_to_wannier();
        //DensityMatrix.go_to_R();
        RK_object.Propagate();

    //Pop.close();
    //Las.close();
    //os_Pos.close();
}

void Simulation::Calculate_Velocity()
{
    for(int ix : {0, 1, 2}){
        Velocity[ix].get_Operator_R() = DensityMatrix.get_Operator_R();
    
        Velocity[ix].get_Operator_R().fill(0.);
        commutator(Velocity[ix].get_Operator_R(), -im, material.r[ix].get_Operator_R(), material.H.get_Operator_R());
        //part with R
        auto& ci = MeshGrid<R>::ConvolutionIndex1[{Velocity[ix].get_Operator_R().get_MeshGrid()->get_id(), 
                                              material.H.get_Operator_R().get_MeshGrid()->get_id(), 
                                              Operator<std::complex<double>>::MeshGrid_Null->get_id()}];
        if(ci.get_Size(0) == 0 ){
            MeshGrid<R>::Calculate_ConvolutionIndex1(*(Velocity[ix].get_Operator_R().get_MeshGrid()), 
                                                     *(material.H.get_Operator_R().get_MeshGrid()), 
                                                     *(Operator<std::complex<double>>::MeshGrid_Null));
        }

        for(int iblock=0; iblock<Velocity[ix].get_Operator_R().get_nblocks(); ++iblock){
            auto& R =(*(Velocity[ix].get_Operator_R().get_MeshGrid()))[iblock].get("Cartesian"); 
            for(int irow=0; irow<Velocity[ix].get_Operator_R().get_nrows(); ++irow){
                for(int icol=0; icol<Velocity[ix].get_Operator_R().get_ncols(); ++icol){
                    if(ci(iblock, 0) != -1){
                        Velocity[ix].get_Operator_R()(iblock, irow, icol) += -im*R[ix]*material.H.get_Operator_R()[ci(iblock,0)](irow, icol);
                    }
                }
            }
        }
    }
}


void Simulation::print_recap()
{
    std::cout << "*************  GRIDS: *************\n";
    std::cout << "Total number of grids: (R) " << MeshGrid<R>::get_counter_id() << " (k) " << MeshGrid<k>::get_counter_id() << std::endl;
    std::cout << "              ";
    std::cout << "| id |";
    std::cout << " size  |\n";
    std::cout << "material grid ";
    std::cout << "|" << std::setw(4) << material.H.get_Operator_R().get_MeshGrid()->get_id() << "|";
    std::cout<<  std::setw(7)  <<  material.H.get_Operator_R().get_MeshGrid()->get_mesh().size() << "|"<< std::endl;
    std::cout << "DM grid (R)   ";
    std::cout << "|" << std::setw(4) << DensityMatrix.get_Operator_R().get_MeshGrid()->get_id() << "|";
    std::cout <<  std::setw(7) << DensityMatrix.get_Operator_R().get_MeshGrid()->get_mesh().size() << "|"<< std::endl;
    std::cout << "DM grid (k)   ";
    std::cout << "|" << std::setw(4) << DensityMatrix.get_Operator_k().get_MeshGrid()->get_id() << "|";
    std::cout <<  std::setw(7) << DensityMatrix.get_Operator_k().get_MeshGrid()->get_mesh().size() << "|"<< std::endl;
    laser.print_info();
}
