#include "Simulation/Simulation.hpp"
#include "core/mpi/Communicator.hpp"

Simulation::Simulation(const std::string& JsonFileName__)
{
#ifdef EDUS_MPI
    if(mpi::Communicator::world().rank() == 0)
#endif
    {
        std::cout << std::setw(125) << "* -> initializing simulation";
        std::cout << "*\n";
    }
    JsonFile = JsonFileName__;
    std::ifstream f(JsonFileName__);

    nlohmann::json data = nlohmann::json::parse(f);

    //read model from tb_file    
    PROFILE("Simulation::Initialize");
    //---------------------------getting info from tb----------------------------------------
    tb_model = data["tb_file"].template get<std::string>();
#ifdef EDUS_MPI
    if(mpi::Communicator::world().rank() == 0)
#endif
    {
        std::cout << std::setw(125) << "* -> initializing material";
        std::cout << "*\n";
    }
    material = Material(tb_model);
    
    //---------------------------------------------------------------------------------------
#ifdef EDUS_MPI
    if(mpi::Communicator::world().rank() == 0)
#endif
    {
        std::cout << std::setw(125) << "* -> initializing lasers";
        std::cout << "*\n";
    }
    for ( int ilaser = 0; ilaser < int( data["lasers"].size() ); ++ilaser ) {
        auto& currentdata = data["lasers"][ilaser];
        Laser laser_;
        laser_.set_InitialTime(0., FemtoSeconds);
        laser_.set_Intensity(currentdata["intensity"][0].template get<double>(), 
                            unit(currentdata["intensity"][1].template get<std::string>()));
        auto freq_wavelength = wavelength_or_frequency(currentdata);
        if(freq_wavelength == "frequency") {
            laser_.set_Omega(currentdata["frequency"][0].template get<double>(), 
                            unit(currentdata["frequency"][1].template get<std::string>()));
        }
        else {
            laser_.set_Lambda(currentdata["wavelength"][0].template get<double>(), 
                            unit(currentdata["wavelength"][1].template get<std::string>()));
        }
        laser_.set_NumberOfCycles(currentdata["cycles"].template get<double>());
        Coordinate pol(currentdata["polarization"][0], currentdata["polarization"][1],
                                          currentdata["polarization"][2]);
        pol = pol/pol.norm(); 
        laser_.set_Polarization(pol);
        setoflaser.push_back(laser_);
    }

    //--------------------------initializing grids and arrays--------------------------------
    auto MasterRgrid = std::make_shared<MeshGrid>(R, data["grid"].template get<std::array<int,3>>());
    FilledBands = data["filledbands"];
    coulomb.set_DoCoulomb(data["coulomb"].template get<bool>());    
    PrintResolution = data["printresolution"].template get<int>();

    InitialTime =    data["initialtime"][0].template get<double>();
    InitialTime = Convert(InitialTime, unit(data["initialtime"][1].template get<std::string>()), AuTime);
    FinalTime =    data["finaltime"][0].template get<double>();
    FinalTime = Convert(FinalTime, unit(data["finaltime"][1].template get<std::string>()), AuTime);
    
    std::array<int, 3> MG_size = {MasterRgrid->get_Size()[0], MasterRgrid->get_Size()[1], MasterRgrid->get_Size()[2]};
    Operator<std::complex<double>>::mpindex.initialize(MG_size);
#ifdef EDUS_MPI
    if(mpi::Communicator::world().rank() == 0)
#endif
    {
        std::cout << std::setw(125) << "* -> initializing DensityMatrix";
        std::cout << "*\n";
    }
    DensityMatrix.initialize_fft(*MasterRgrid, material.H.get_Operator_R().get_nrows());
    Operator<std::complex<double>>::SpaceOfPropagation = SpaceOfPropagation;
#ifdef EDUS_MPI
    if(mpi::Communicator::world().rank() == 0)
#endif
    {
        std::cout << std::setw(125) << "* -> H to k";
        std::cout << "*\n";
    }

    material.H.dft(DensityMatrix.get_FT_meshgrid_k().get_mesh(), +1);
    for(auto ix : {0,1,2}) {
        material.r[ix].dft(DensityMatrix.get_FT_meshgrid_k().get_mesh(), +1);
        material.r[ix].get_Operator(Space::k).make_hermitian();
    }
#ifdef EDUS_MPI
    if(mpi::Communicator::world().rank() == 0)
#endif
    {
        std::cout << std::setw(125) << "* -> initialize fft";
        std::cout << "*\n";
    }
    H.initialize_fft(DensityMatrix);
#ifdef EDUS_MPI
    if(mpi::Communicator::world().rank() == 0)
#endif
    {
        std::cout << std::setw(125) << "* -> eigensolver";
        std::cout << "*\n";
    }
    SettingUp_EigenSystem();
    auto& Uk = Operator<std::complex<double>>::EigenVectors;

    //---------------------------------------------------------------------------------------

    //------------------------setting up TD equations----------------------------------------
    #include "Functional_InitialCondition.hpp"
    #include "Functional_SourceTerm.hpp"
    RK_object.initialize(DensityMatrix, 
                        InitialCondition, SourceTerm);
    RK_object.set_InitialTime(InitialTime);
    RK_object.set_ResolutionTime( Convert(data["dt"][0].template get<double>(), 
                                          unit(data["dt"][1].template get<std::string>()), 
                                        AuTime ));

    kgradient.initialize(*(DensityMatrix.get_Operator(R).get_MeshGrid()));
    coulomb.initialize(material.H.get_Operator_R().get_nrows(), DensityMatrix.get_Operator(R).get_MeshGrid(), material.r);


    print_recap();
    //---------------------------------------------------------------------------------------

    //print DM in R to prove it decays and is zero for large R
    std::stringstream rank;
#ifdef EDUS_MPI
    rank << "DM" << mpi::Communicator::world().rank() << ".txt";
#else
    rank << "DM0.txt";
#endif
    DensityMatrix.go_to_R();
    std::ofstream os;
    os.open(rank.str());
    auto Rgamma_centered = get_GammaCentered_grid(*DensityMatrix.get_Operator_R().get_MeshGrid());
    for(int iR_loc=0; iR_loc< DensityMatrix.get_Operator_R().get_nblocks(); ++iR_loc){
        //os << DensityMatrix.get_Operator_R()[i] << std::endl;
        auto iR_glob = DensityMatrix.mpindex.loc1D_to_glob1D(iR_loc);
        os << Rgamma_centered[iR_glob].norm();
        os << " " << std::abs(max(DensityMatrix.get_Operator_R()[iR_loc])) << std::endl;
    }
    os.close();

    coulomb.set_DM0(DensityMatrix);

    Calculate_Velocity();
    DensityMatrix.go_to_k();
    assert(DensityMatrix.get_Operator_k().is_hermitian());
    //---------------------------------------------------------------------------------------

    //---------------------open files---------------------------------------------------------
    os_Pop.open("Population.txt");
    os_Laser.open("Laser.txt");
    os_VectorPot.open("Laser_A.txt");
    os_Time.open("Time.txt");
    os_Velocity.open("Velocity.txt");
    //---------------------------------------------------------------------------------------
}




bool Simulation::PrintObservables(const double& time) const
{
    return ( int( round( time/RK_object.get_ResolutionTime() ) ) % PrintResolution == 0 );
}

void Simulation::SettingUp_EigenSystem()
{
    PROFILE("Simulation::SetEigensystem");

    //------------------------go to H(k)--------------------------------------------
    auto& MasterkGrid = DensityMatrix.get_Operator_k().get_MeshGrid()->get_mesh();
    //material.H.dft(MasterkGrid, +1);
    //------------------------------------------------------------------------------

    //--------------------solve eigen problem---------------------------------------
    auto& Uk = Operator<std::complex<double>>::EigenVectors;
    auto& UkDagger = Operator<std::complex<double>>::EigenVectors_dagger;

    material.H.get_Operator_k().diagonalize(Band_energies, Uk);
    //------------------------------------------------------------------------------

    //-------------------get U dagger-----------------------------------------------
    UkDagger.initialize(k, Uk.get_nblocks(), Uk.get_ncols(), Uk.get_nrows());
    for(int ik=0; ik<UkDagger.get_nblocks(); ++ik){
        for(int ir=0; ir<UkDagger.get_nrows(); ++ir){
            for(int ic=0; ic<UkDagger.get_ncols(); ++ic){
                UkDagger[ik](ir, ic) = conj(Uk[ik](ic,ir));
            }
        }
    }
    //------------------------------------------------------------------------------
};   


void Simulation::Calculate_TDHamiltonian(const double& time, const bool& erase_H)
{
#ifdef EDUS_TIMERS
    PROFILE("SourceTerm::Calculate_TDHamiltonian");
#endif
    H.go_to_k();
    //--------------------get aliases for nested variables--------------------------------
    auto& H_ = H.get_Operator(SpaceOfPropagation);
    auto& H0_ = material.H.get_Operator(SpaceOfPropagation);
    auto& x_ = material.r[0].get_Operator(SpaceOfPropagation);
    auto& y_ = material.r[1].get_Operator(SpaceOfPropagation);
    auto& z_ = material.r[2].get_Operator(SpaceOfPropagation);
    auto las = setoflaser(time).get("Cartesian");
    
    //auto& ci = MeshGrid::ConvolutionIndex[{H0_.get_MeshGrid()->get_id(), 
    //                                          H_.get_MeshGrid()->get_id(), 
    //                                          Operator<std::complex<double>>::MeshGrid_Null->get_id()}];
    //-----------------------------------------------------------------------------------

    //--------------------------do initializations---------------------------------------
    if(erase_H) {
        H_.fill(0.);
    }

    //if(ci.get_Size(0) == 0 && SpaceOfPropagation == R) {
    //    MeshGrid::Calculate_ConvolutionIndex( *(H0_.get_MeshGrid()), 
    //                                              *(H_.get_MeshGrid()), 
    //                                              *(Operator<std::complex<double>>::MeshGrid_Null));
    //}
    //---------------------------------------------------------------------------------

    //------------------------H(R) = H0(R) + E.r(R)-------------------------------------
    #pragma omp parallel for schedule(static) collapse(3)
    for(int iblock=0; iblock<H0_.get_nblocks(); ++iblock){
        for(int irow=0; irow<H0_.get_nrows(); ++irow){
            for(int icol=0; icol<H0_.get_ncols(); ++icol){
                //auto Hblock = ( SpaceOfPropagation == k ? iblock : ci(iblock, 0) );
                //H_(Hblock, irow, icol) = H0_(iblock, irow, icol)
                H_(iblock, irow, icol) += H0_(iblock, irow, icol)
                                        + las[0]*x_(iblock, irow, icol)
                                        + las[1]*y_(iblock, irow, icol)
                                        + las[2]*z_(iblock, irow, icol);
            }
        }
    }
    //---------------------------------------------------------------------------------
}

void Simulation::Propagate()
{
    PROFILE("Simulation::Propagate");
    int iFinalTime = int((FinalTime - InitialTime)/RK_object.get_ResolutionTime())+2;
#ifdef EDUS_MPI
    if(mpi::Communicator::world().rank() == 0)
#endif
    {
        std::cout << "*******************************************************    PROPAGATION     ***************************************************\n";
        std::cout << "*   Initial time:       *     ";
        std::cout << std::left << std::scientific << std::setw(10) << std::setprecision(4) << InitialTime;
        std::cout << std::left << std::setw(15) << " a.u.";
        std::cout << std::left << std::scientific << std::setw(10) << std::setprecision(4) << Convert(InitialTime, AuTime, FemtoSeconds);
        std::cout << std::left << std::setw(60) << " fs" << "*\n";
        std::cout << "*   Final time:         *     ";
        std::cout << std::left << std::scientific << std::setw(10) << std::setprecision(4) << FinalTime;
        std::cout << std::left << std::setw(15) << " a.u.";
        std::cout << std::left << std::scientific << std::setw(10) << std::setprecision(4) << Convert(FinalTime, AuTime, FemtoSeconds);
        std::cout << std::left << std::setw(60) << " fs" << "*\n";
        std::cout << "*   Number of steps:    *     ";
        std::cout << std::left << std::setw(95) << iFinalTime <<  "*\n";
    }
    for( int it = 0; it < iFinalTime; ++it ) {
        if(it%100 == 0) {
#ifdef EDUS_MPI
    if(mpi::Communicator::world().rank() == 0)
#endif
            {
            std::cout << "*   it:                 *     ";
            std::cout << std::setw(5) << std::left <<  it;
            std::cout << " / "; 
            std::cout << std::setw(10) << std::left << iFinalTime;
            std::cout << std::fixed << std::left << std::setw(8) << std::setprecision(4) << 100*double(it)/iFinalTime << "%";
            std::cout << std::right << std::setw(69) << "*" << std::endl;
            }
        }
        do_onestep();
    }
#ifdef EDUS_MPI
    if(mpi::Communicator::world().rank() == 0)
#endif
    std::cout << "******************************************************************************************************************************\n";
}

void Simulation::do_onestep()
{
    auto CurrentTime = RK_object.get_CurrentTime();
    //------------------------Print population-------------------------------------
    
    if( PrintObservables( CurrentTime) ){
        //print time 
        os_Time << CurrentTime << std::endl;
        //print laser
        os_Laser << setoflaser(RK_object.get_CurrentTime()).get("Cartesian");
        os_VectorPot << setoflaser.VectorPotential(RK_object.get_CurrentTime()).get("Cartesian");

        Print_Population();
        Print_Velocity();
    }
    //------------------------------------------------------------------------------
    RK_object.Propagate();
}


void Simulation::Print_Population()
{
    static Operator<std::complex<double>> aux_DM;
    aux_DM = DensityMatrix;
    aux_DM.go_to_bloch();

    auto Population = TraceK(aux_DM.get_Operator(Space::k));
#ifdef EDUS_MPI
    if( kpool_comm.rank() == 0 )
#endif
    {
        for(int ibnd=0; ibnd < int( Population.size() ); ibnd++) {
            Population[ ibnd ] /= double(aux_DM.get_Operator_k().get_MeshGrid()->get_TotalSize());        
        }
        for(int ibnd = 0; ibnd < int( Population.size() ); ibnd++){
            if( ibnd < FilledBands ) Population[ibnd] = 1.-Population[ibnd];
            os_Pop << std::setw(20) << std::setprecision(6) << Population[ibnd].real();
            os_Pop << " ";
        }
        os_Pop << std::endl;
    }
}


void Simulation::Calculate_Velocity()
{
    Calculate_TDHamiltonian(-6000, true);
    H.go_to_R();
    std::vector<Coordinate> direction(3);
    direction[0].initialize(1,0,0);
    direction[1].initialize(0,1,0);
    direction[2].initialize(0,0,1);

    for(int ix : {0, 1, 2}){
        Velocity[ix].initialize_fft(DensityMatrix);
        Velocity[ix].lock_space(k);
        Velocity[ix].get_Operator_k().fill(0.);
        commutator(Velocity[ix].get_Operator_k(), -im, material.r[ix].get_Operator_k(), H.get_Operator_k());
        Velocity[ix].go_to_R();
        //part with R
        kgradient.Calculate(1., Velocity[ix].get_Operator_R(), H.get_Operator_R(), direction[ix], false);
        Velocity[ix].go_to_k();
    }
}


void Simulation::Print_Velocity()
{
    std::array<std::complex<double>, 3> v = {0., 0., 0.};
    
    auto& DMK = DensityMatrix.get_Operator_k();
    static BlockMatrix<std::complex<double>> temp(k, DMK.get_nblocks(), DMK.get_nrows(), DMK.get_ncols());

    for(auto ix : {0, 1, 2}){
        temp.fill(0.);
        multiply(temp, 1.+im*0., Velocity[ix].get_Operator_k(), DMK );
        auto ToSum = TraceK(temp);
        for(int ibnd=0; ibnd<int(ToSum.size()); ++ibnd) {
            v[ix] += ToSum[ibnd];
        }
        v[ix] /= DMK.get_MeshGrid()->get_TotalSize();
    }
#ifdef EDUS_MPI
    if( kpool_comm.rank() == 0 )
#endif
    {
        os_Velocity << std::setw(20) << std::setprecision(8) << v[0].real();
        os_Velocity << std::setw(20) << std::setprecision(8) << v[0].imag();
        os_Velocity << std::setw(20) << std::setprecision(8) << v[1].real();
        os_Velocity << std::setw(20) << std::setprecision(8) << v[1].imag();
        os_Velocity << std::setw(20) << std::setprecision(8) << v[2].real();
        os_Velocity << std::setw(20) << std::setprecision(8) << v[2].imag();
        os_Velocity << std::endl;
    }
}

void Simulation::print_recap()
{
    //std::cout << "*************  GRIDS: *************\n";
    //std::cout << "Total number of grids:  " << MeshGrid::get_counter_id() << std::endl;
    //std::cout << "              ";
    //std::cout << "| id |";
    //std::cout << " size  |\n";
    //std::cout << "material grid ";
    //std::cout << "|" << std::setw(4) << material.H.get_Operator_R().get_MeshGrid()->get_id() << "|";
    //std::cout<<  std::setw(7)  <<  material.H.get_Operator_R().get_MeshGrid()->get_mesh().size() << "|"<< std::endl;
    //std::cout << "DM grid (R)   ";
    //std::cout << "|" << std::setw(4) << DensityMatrix.get_Operator_R().get_MeshGrid()->get_id() << "|";
    //std::cout <<  std::setw(7) << DensityMatrix.get_Operator_R().get_MeshGrid()->get_mesh().size() << "|"<< std::endl;
    //std::cout << "DM grid (k)   ";
    //std::cout << "|" << std::setw(4) << DensityMatrix.get_Operator_k().get_MeshGrid()->get_id() << "|";
    //std::cout <<  std::setw(7) << DensityMatrix.get_Operator_k().get_MeshGrid()->get_mesh().size() << "|"<< std::endl;
    //"tb_file": "/home/gcistaro/EDUS/tb_models/hBN_gap7.25eV_a2.5A",
    //"dt": [0.1, "autime"],
    //"solver": "RungeKutta",
    //"printresolution": 6,
    //"coulomb": false,
#ifdef EDUS_MPI
    if(mpi::Communicator::world().rank() == 0)
#endif
    {
        std::cout << "*******************************************************    INPUT RECAP     ***************************************************\n";
        std::cout << "*   input file:         *     ";
        std::cout << std::left << std::setw(95) << JsonFile << "*\n";
        std::cout << "*   tb_model  :         *     ";
        std::cout << std::left << std::setw(95) << tb_model << "*\n";
        std::cout << "*   grid:               *     [";
        std::cout << std::right << std::setw(4) << DensityMatrix.get_Operator_R().get_MeshGrid()->get_Size()[0];
        std::cout << ", ";
        std::cout << std::right << std::setw(4) << DensityMatrix.get_Operator_R().get_MeshGrid()->get_Size()[1];
        std::cout << ", ";
        std::cout << std::right << std::setw(4) << DensityMatrix.get_Operator_R().get_MeshGrid()->get_Size()[2];
        std::cout << std::left << std::setw(78) << "]" << "*\n";     
        std::cout << "*   Resolution time:    *     ";
        std::cout << std::left << std::scientific << std::setw(10) << std::setprecision(4) << RK_object.get_ResolutionTime();
        std::cout << std::left << std::setw(15) << " a.u.";
        std::cout << std::left << std::scientific << std::setw(10) << std::setprecision(4) << Convert(RK_object.get_ResolutionTime(), AuTime, FemtoSeconds);
        std::cout << std::left << std::setw(60) << " fs" << "*\n";
        std::cout << "*   Print Resolution:   *     ";
        std::cout << std::left << std::setw(95) << PrintResolution <<  "*\n";
        std::cout << "*   Coulomb:            *     ";
        std::cout << std::left << std::setw(95)  << (coulomb.get_DoCoulomb() ? "True" : "False") << "*\n";
        std::cout << "******************************************************************************************************************************\n";
        std::cout << "*******************************************************   WANNIER     ********************************************************\n";
        std::cout << "*   #R points :         *     ";
        std::cout << std::setw(95) << std::left<< material.H.get_Operator_R().get_nblocks();
        std::cout << "*\n";
        auto& A = Coordinate::get_Basis(LatticeVectors(R)).get_M();
        std::cout << "*             :         *   ";
        std::cout << std::scientific << std::right << std::setw(14) << std::setprecision(6) <<  A(0,0);
        std::cout << std::scientific << std::right << std::setw(14) << std::setprecision(6) <<  A(0,1);
        std::cout << std::scientific << std::right << std::setw(14) << std::setprecision(6) <<  A(0,2);
        std::cout << std::setw(56) << std::right << "*" << std::endl;
        std::cout << "*      A      :         *   ";
        std::cout << std::scientific << std::right << std::setw(14) << std::setprecision(6) <<  A(1,0);
        std::cout << std::scientific << std::right << std::setw(14) << std::setprecision(6) <<  A(1,1);
        std::cout << std::scientific << std::right << std::setw(14) << std::setprecision(6) <<  A(1,2);
        std::cout << std::setw(56) << std::right << "*" << std::endl;
        std::cout << "*             :         *   ";
        std::cout << std::scientific << std::right << std::setw(14) << std::setprecision(6) <<  A(2,0);
        std::cout << std::scientific << std::right << std::setw(14) << std::setprecision(6) <<  A(2,1);
        std::cout << std::scientific << std::right << std::setw(14) << std::setprecision(6) <<  A(2,2);
        std::cout << std::setw(56) << std::right << "*" << std::endl;
        std::cout << "*";
        std::cout << std::setw(124) << " ";
        std::cout << "*\n";
        auto& B = Coordinate::get_Basis(LatticeVectors(k)).get_M();
        std::cout << "*             :         *   ";
        std::cout << std::scientific << std::right << std::setw(14) << std::setprecision(6) <<  B(0,0);
        std::cout << std::scientific << std::right << std::setw(14) << std::setprecision(6) <<  B(0,1);
        std::cout << std::scientific << std::right << std::setw(14) << std::setprecision(6) <<  B(0,2);
        std::cout << std::setw(56) << std::right << "*" << std::endl;
        std::cout << "*      B      :         *   ";
        std::cout << std::scientific << std::right << std::setw(14) << std::setprecision(6) <<  B(1,0);
        std::cout << std::scientific << std::right << std::setw(14) << std::setprecision(6) <<  B(1,1);
        std::cout << std::scientific << std::right << std::setw(14) << std::setprecision(6) <<  B(1,2);
        std::cout << std::setw(56) << std::right << "*" << std::endl;
        std::cout << "*             :         *   ";
        std::cout << std::scientific << std::right << std::setw(14) << std::setprecision(6) <<  B(2,0);
        std::cout << std::scientific << std::right << std::setw(14) << std::setprecision(6) <<  B(2,1);
        std::cout << std::scientific << std::right << std::setw(14) << std::setprecision(6) <<  B(2,2);
        std::cout << std::setw(56) << std::right << "*" << std::endl;
        std::cout << "******************************************************************************************************************************\n";
        for( int ilaser=0; ilaser<int(setoflaser.size()); ++ilaser) {
            setoflaser[ilaser].print_info();
        }
    }
}


int Simulation::get_it(const double& time) const
{
    return int(time/RK_object.get_ResolutionTime());
}




void Simulation::print_grids()
{
    auto& MasterRgrid = DensityMatrix.get_Operator_R().get_MeshGrid();
    auto& kgrid = DensityMatrix.get_Operator_k().get_MeshGrid();
    auto Rfft = std::make_shared<MeshGrid>(fftPair(*kgrid));//Rspace as fourier space of kgrid
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
    auto MeshRfft = fftPair(*kgrid);
    os.open("R.txt");
    for (auto &R : Rfft->get_mesh())
    {
        os << R.get("Cartesian");
    }
    os.close();
}




std::string wavelength_or_frequency(const nlohmann::json& data)
{
    bool is_frequency = ( data.find("frequency") != data.end() );
    if( is_frequency) return "frequency";
    bool is_wavelength = ( data.find("wavelength") != data.end() );
    if( !is_wavelength ) {
        std::cout << "You must specify frequency *xor* wavelength!\n";
        exit(1);
    }
    return "wavelength";
}