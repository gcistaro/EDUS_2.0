#include "Simulation/Simulation.hpp"

Simulation::Simulation(const std::string& FileName, const double& Radius)
{
    PROFILE_START("Simulation::Initialize");
    material = Material(FileName);
    laser.set_Intensity(0., Wcm2);//1.e+16, Wcm2);
    laser.set_InitialTime(0., FemtoSeconds);
    laser.set_Lambda(800, NanoMeters);
    laser.set_NumberOfCycles(5);
    laser.set_Polarization(Coordinate<k>(1,1,1));
    laser.print_info();
    //Preparing grids and calculating eigenvectors
    //std::array<int,3> Size = {20,20,1};
    //auto MasterRgrid = std::make_shared<MeshGrid<R>>(Size);//MeshGrid<R>(Radius));//spherical grid for exponentially decaying DM
    auto MasterRgrid = std::make_shared<MeshGrid<R>>(Radius);//spherical grid for exponentially decaying DM
    //get U and Udagger from Hamiltonian
    DensityMatrix.initialize_fft(*MasterRgrid, material.H.get_Operator_R().get_nrows());

    H = DensityMatrix; //to initialize dimensions
    std::cout << std::endl;
    SettingUp_EigenSystem();
    auto& Uk = Operator<std::complex<double>>::EigenVectors;
    auto& UkDagger = Operator<std::complex<double>>::EigenVectors_dagger;
   
    #include "Functionals.hpp"
    //RK_object.initialize(DensityMatrix.get_Operator_k(), 
    //                    InitialCondition_k, SourceTerm_k);
    RK_object.initialize(DensityMatrix.get_Operator_R(), 
                        InitialCondition, SourceTerm);

    //print DM in R 
    std::ofstream os;
    os.open("DM.txt");
    for(int i=0; i< DensityMatrix.get_Operator_R().get_MeshGrid()->get_mesh().size(); i++){
        os << DensityMatrix.get_Operator_R().get_MeshGrid()->get_mesh()[i].norm() << " " << std::abs(max(DensityMatrix.get_Operator_R()[i])) << std::endl;
    }
    os.close();

    RK_object.set_InitialTime(0.);
    RK_object.set_ResolutionTime(0.01);
    print_recap();
    PROFILE_STOP("Simulation::Initialize");
    this->Propagate();
}

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
    HR.fill(0);
    //TODO: CREATE A COMMONG MESHGRID_NULL!!!! WASTE!!!!!!!!!!!!
    mdarray<double,2> bare_mg({1,3});
    bare_mg.fill(0);
    auto MeshGrid_Null = std::make_shared<MeshGrid<R>>(bare_mg, "Cartesian");
    ////////////////////////////////////////////////////////////
    auto& ci = MeshGrid<R>::ConvolutionIndex1[{HR.get_MeshGrid()->get_id(), 
                                              H0R.get_MeshGrid()->get_id(), 
                                              MeshGrid_Null->get_id()}];
    if(ci.get_Size(0) == 0 ){
        MeshGrid<R>::Calculate_ConvolutionIndex1(*(HR.get_MeshGrid()) , *(H0R.get_MeshGrid()), *MeshGrid_Null);
    }

    //non-easy part
    for(int iblock=0; iblock<HR.get_nblocks(); ++iblock){
        for(int irow=0; irow<H.get_Operator_R().get_nrows(); ++irow){
            for(int icol=0; icol<H.get_Operator_R().get_ncols(); ++icol){
                HR(iblock, irow, icol) = H0R(ci(iblock, 0), irow, icol);
                                        + las[0]*xR(ci(iblock, 0), irow, icol);
                                        + las[1]*yR(ci(iblock, 0), irow, icol);
                                        + las[2]*zR(ci(iblock, 0), irow, icol);
            }
        }
    }
    //R operator acts as R-R'rho(R')
    for(int iblock=0; iblock<H0R.get_nblocks(); ++iblock){
        auto& R = HR.get_MeshGrid()->get_mesh()[iblock].get("Cartesian");
        for(int irow=0; irow<H.get_Operator_R().get_nrows(); ++irow){
            HR(iblock, irow, irow) += las[0]*R[0];
            HR(iblock, irow, irow) += las[1]*R[1];
            HR(iblock, irow, irow) += las[2]*R[2];
        }
    }
}

void Simulation::Propagate()
{
    for(int i=0; i<10; i++){
        //std::cout << "i " << i << std::endl;
        PROFILE("RK_Propagate");
        DensityMatrix.go_to_k();
        DensityMatrix.go_to_bloch();
        
        std::vector<std::complex<double>> Population(DensityMatrix.get_Operator_k().get_nrows(), 0.);
        for(int ik=0; ik<DensityMatrix.get_Operator_k().get_nblocks(); ik++){
            for(int ibnd=0; ibnd<Population.size(); ibnd++){
                Population[ibnd] += DensityMatrix.get_Operator_k()[ik](ibnd,ibnd);
            }
        }
        std::cout << std::endl << std::endl << std::endl << "POPULATION: ";
        for(int ibnd=0; ibnd<Population.size(); ibnd++){
            std::cout << Population[ibnd] << ' ';
        }

        //std::cout << "Density Matrix R : " << DensityMatrix.get_Operator_R() << std::endl;
        //std::cout << std::endl;
        //exit(0);
        //std::cout << i << std::endl;
        DensityMatrix.go_to_wannier();
        DensityMatrix.go_to_R();
        RK_object.Propagate();
        std::cout << std::endl << std::endl;
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

