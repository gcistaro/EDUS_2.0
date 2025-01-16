#include <filesystem>
#include "Simulation/Simulation.hpp"
#include "core/mpi/Communicator.hpp"

Simulation::Simulation(const std::string& JsonFileName__)
{
    output::print("-> initializing simulation");

    JsonFile = JsonFileName__;
    std::ifstream f(JsonFileName__);

    nlohmann::json data = nlohmann::json::parse(f);

    PROFILE("Simulation::Initialize");

    //---------------------------getting info from tb----------------------------------------
    tb_model = data["tb_file"].template get<std::string>();
    output::print("-> initializing material");
    material = Material(tb_model);
    
    //---------------------------------------------------------------------------------------
    output::print("-> initializing lasers");

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
            currentdata["wavelength"] = {laser_.get_Lambda()};
        }
        else {
            laser_.set_Lambda(currentdata["wavelength"][0].template get<double>(), 
                            unit(currentdata["wavelength"][1].template get<std::string>()));
            currentdata["frequency"] = {laser_.get_Omega()};
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
    Operator<std::complex<double>>::mpindex.initialize(MG_size, material.H.get_Operator_R().get_nrows()*material.H.get_Operator_R().get_nrows());

    output::print("-> initializing DensityMatrix");
    DensityMatrix.initialize_fft(*MasterRgrid, material.H.get_Operator_R().get_nrows());
    Operator<std::complex<double>>::SpaceOfPropagation = SpaceOfPropagation;
    output::print("-> H to k");
    material.H.dft(DensityMatrix.get_FT_meshgrid_k().get_mesh(), +1);
    for(auto ix : {0,1,2}) {
        material.r[ix].dft(DensityMatrix.get_FT_meshgrid_k().get_mesh(), +1);
        material.r[ix].get_Operator(Space::k).make_hermitian();
    }

    output::print("-> initialize fft");
    H.initialize_fft(DensityMatrix);

    output::print("-> solve eigensystem");
    SettingUp_EigenSystem();
    auto& Uk = Operator<std::complex<double>>::EigenVectors;

    //---------------------------------------------------------------------------------------

    //------------------------setting up TD equations----------------------------------------
    #include "Functional_InitialCondition.hpp"
    #include "Functional_SourceTerm.hpp"
    DEsolver_DM.initialize(DensityMatrix, 
                        InitialCondition, SourceTerm, 
                        solver.at(data["solver"].template get<std::string>()),
                        data["order"].template get<int>()  );
    DEsolver_DM.set_InitialTime(InitialTime);
    DEsolver_DM.set_ResolutionTime( Convert(data["dt"][0].template get<double>(), 
                                          unit(data["dt"][1].template get<std::string>()), 
                                        AuTime ));

    kgradient.initialize(*(DensityMatrix.get_Operator(R).get_MeshGrid()));
    coulomb.initialize(material.H.get_Operator_R().get_nrows(), DensityMatrix.get_Operator(R).get_MeshGrid(), material.r);
    DensityMatrix.go_to_R();
    coulomb.set_DM0(DensityMatrix);


    print_recap();
    //---------------------------------------------------------------------------------------

    //print DM in R to prove it decays and is zero for large R
    std::filesystem::create_directories(std::string(std::filesystem::current_path())+"/Output");
    std::filesystem::current_path(std::string(std::filesystem::current_path())+"/Output");


    std::stringstream rank;
#ifdef EDUS_MPI
    rank << "DM" << mpi::Communicator::world().rank() << ".txt";
#else
    rank << "DM0.txt";
#endif
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

    // ---------------------- create recap file ---------------------------------------------
    data["num_bands"] = {material.H.get_Operator_R().get_nrows()};
    data["num_kpoints"] = {material.H.get_Operator_k().get_MeshGrid()->get_TotalSize()};

    
    mdarray<double,2> bare_k({material.H.get_Operator_k().get_MeshGrid()->get_TotalSize(),3});
    for( int ik=0; ik < material.H.get_Operator_k().get_MeshGrid()->get_mesh().size(); ++ik ) {
        for(auto& ix : {0, 1, 2}) {
            bare_k(ik,ix) = (*(material.H.get_Operator_k().get_MeshGrid()))[ik].get(LatticeVectors(k))[ix];
        }
    }
    data["kpoints"] = bare_k;
    data["A"] = Coordinate::get_Basis(LatticeVectors(R)).get_M();
    data["B"] = Coordinate::get_Basis(LatticeVectors(k)).get_M();

#ifdef EDUS_HDF5
    std::string name = "output.h5";
    if (!file_exists(name)) {
        HDF5_tree(name, hdf5_access_t::truncate);
    }
    HDF5_tree fout(name, hdf5_access_t::truncate);
    dump_json_in_h5( data, name );
    mpi::Communicator::world().barrier();
#endif

    //---------------------------------------------------------------------------------------
}




bool Simulation::PrintObservables(const double& time) const
{
    return ( int( round( time/DEsolver_DM.get_ResolutionTime() ) ) % PrintResolution == 0 );
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
    // == /* H_ = H0_ + las_x \cdot x */
    // == SumWithProduct(H_, 1., H0_, las[0], x_);
    // == /* H_ = H_ + las_y \cdot y */
    // == SumWithProduct(H_, 1., H_, las[1], y_);
    // == /* H_ = H_ + las_z \cdot z */
    // == SumWithProduct(H_, 1., H0_, las[1], y_);

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
    int iFinalTime = int((FinalTime - InitialTime)/DEsolver_DM.get_ResolutionTime())+2;

    /* print info on output */
    output::title("PROPAGATION");
    output::print("Initial time:       *", InitialTime, " a.u.", Convert(InitialTime, AuTime, FemtoSeconds), " fs");
    output::print("Final time:         *", FinalTime, " a.u.", Convert(FinalTime, AuTime, FemtoSeconds), " fs");
    output::print("Number of steps:    *", iFinalTime);


    /* do steps */
    for( int it = 0; it < iFinalTime; ++it ) {
        if(it%100 == 0) {
            output::print( "it:                 *", it, " / ", iFinalTime, 100*double(it)/iFinalTime, " %");
        }
        do_onestep();
    }

    output::print(std::string(output::linesize-5, '*'));
}

void Simulation::do_onestep()
{
    auto CurrentTime = DEsolver_DM.get_CurrentTime();
    //------------------------Print population-------------------------------------
    
    if( PrintObservables( CurrentTime ) ){
#ifdef EDUS_HDF5
        std::string name_ = "output.h5";
        auto node = get_it_sparse(CurrentTime);

        HDF5_tree fout(name_, hdf5_access_t::truncate);

        fout.create_node(node);
        fout.write("time_au", CurrentTime);
        Calculate_TDHamiltonian(CurrentTime, true);
        H.get_Operator_k().write_h5(name_, node, nodename::fullH);
        DensityMatrix.go_to_R(true);
        H.go_to_R(false);
        coulomb.EffectiveHamiltonian( H, DensityMatrix, true); 
        H.go_to_k(true);
        H.get_Operator_k().write_h5(name_, node, "SelfEnergy");
        DensityMatrix.get_Operator_k().write_h5(name_, node, nodename::DMk);
        DensityMatrix.go_to_k(false);
#endif 

        //print time 
        os_Time << CurrentTime << std::endl;
        //print laser
        os_Laser << setoflaser(DEsolver_DM.get_CurrentTime()).get("Cartesian");
        os_VectorPot << setoflaser.VectorPotential(DEsolver_DM.get_CurrentTime()).get("Cartesian");

        Print_Population();
        Print_Velocity();
    }
    //------------------------------------------------------------------------------
    DEsolver_DM.Propagate();
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
    std::stringstream title;
    title << std::string(5, ' ') <<  "INPUT RECAP" << std::string(5, ' ');
    int num_stars = (output::linesize - title.str().length())/2-2;

    output::print(std::string(num_stars,'*'), title.str(), std::string(num_stars,'*'));
    output::print("input file:         *", std::string(8, ' '), JsonFile );
    output::print("tb_model  :         *", std::string(8, ' '), tb_model);
    output::print("grid      :         *", std::string(8, ' '), "[", 
                                        DensityMatrix.get_Operator_R().get_MeshGrid()->get_Size()[0], ", ",
                                        DensityMatrix.get_Operator_R().get_MeshGrid()->get_Size()[1], ", ",
                                        DensityMatrix.get_Operator_R().get_MeshGrid()->get_Size()[2], "]");
    output::print("Resolution time:    *", DEsolver_DM.get_ResolutionTime(), " a.u.", Convert(DEsolver_DM.get_ResolutionTime(), AuTime, FemtoSeconds), " fs");
    output::print("PrintResolution:    *", PrintResolution);
    output::print("Coulomb        :    *", std::string(8, ' '), (coulomb.get_DoCoulomb() ? "True" : "False"));
    output::stars();

    output::stars();
    title.clear();
    title << std::string(5, ' ') <<  "WANNIER" << std::string(5, ' ');
    num_stars = (output::linesize - title.str().length())/2-2;
    output::print(std::string(num_stars,'*'), title.str(), std::string(num_stars,'*'));

    output::print("#R points :         *", material.H.get_Operator_R().get_nblocks());
    auto& A = Coordinate::get_Basis(LatticeVectors(R)).get_M();
    output::print("          :         *", A(0,0), A(0,1), A(0,2));
    output::print("    A     :         *", A(1,0), A(1,1), A(1,2));
    output::print("          :         *", A(2,0), A(2,1), A(2,2));
    auto& B = Coordinate::get_Basis(LatticeVectors(k)).get_M();
    output::print("          :         *", B(0,0), B(0,1), B(0,2));
    output::print("    B     :         *", B(1,0), B(1,1), B(1,2));
    output::print("          :         *", B(2,0), B(2,1), B(2,2));
    output::stars();

    for( int ilaser=0; ilaser<int(setoflaser.size()); ++ilaser) {
        setoflaser[ilaser].print_info();
    }
}



int Simulation::get_it(const double& time) const
{
    return int(time/DEsolver_DM.get_ResolutionTime());
}

int Simulation::get_it_sparse(const double& time) const
{
    return int(round(time/DEsolver_DM.get_ResolutionTime()/PrintResolution));
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
        throw std::runtime_error("You must specify frequency *xor* wavelength!");
    }
    return "wavelength";
}
