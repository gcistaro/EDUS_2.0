/*

template<class T>
Simulation::Simulation(const std::string& FileName, const T& arg_meshinit)
{
    PROFILE("Simulation::Initialize");
    //---------------------------getting info from tb----------------------------------------
    material = Material(FileName);
    //---------------------------------------------------------------------------------------
    
    //------------------------------initializing laser---------------------------------------
    laser.set_InitialTime(0., FemtoSeconds);
    //laser.set_Lambda(3000, NanoMeters);
    laser.set_Intensity(1.e+05, Wcm2);
    laser.set_Lambda(232.368, NanoMeters);
    laser.set_NumberOfCycles(1);
    laser.set_Polarization(Coordinate(1,0,0));
    //---------------------------------------------------------------------------------------

    //--------------------------initializing grids and arrays--------------------------------
    auto MasterRgrid = std::make_shared<MeshGrid>(R, arg_meshinit);
    std::array<int, 3> MG_size = {MasterRgrid->get_Size()[0], MasterRgrid->get_Size()[1], MasterRgrid->get_Size()[2]};
    Operator<std::complex<double>>::mpindex.initialize(MG_size);
    DensityMatrix.initialize_fft(*MasterRgrid, material.H.get_Operator_R().get_nrows());
    Operator<std::complex<double>>::SpaceOfPropagation = SpaceOfPropagation;

    material.H.dft(DensityMatrix.get_FT_meshgrid_k().get_mesh(), +1);
    material.r[0].dft(DensityMatrix.get_FT_meshgrid_k().get_mesh(), +1);
    material.r[1].dft(DensityMatrix.get_FT_meshgrid_k().get_mesh(), +1);
    material.r[2].dft(DensityMatrix.get_FT_meshgrid_k().get_mesh(), +1);
    std::stringstream rank_H;
    rank_H << "H" << mpi::Communicator::world().rank() << ".txt";
    std::ofstream os_H;
    os_H.open(rank_H.str());
    for(int i=0; i< H.get_Operator_k().get_nblocks(); ++i){
        os_H << material.H.get_Operator_k()[i] << std::endl;
    }
    std::stringstream rank_x;
    rank_x << "x" << mpi::Communicator::world().rank() << ".txt";
    std::ofstream os_x;
    os_x.open(rank_x.str());
    for(int i=0; i< material.r[0].get_Operator_k().get_nblocks(); ++i){
        os_x << material.r[0].get_Operator_k()[i] << std::endl;
    }
    std::stringstream rank_y;
    rank_y << "y" << mpi::Communicator::world().rank() << ".txt";
    std::ofstream os_y;
    os_y.open(rank_y.str());
    for(int i=0; i< material.r[1].get_Operator_k().get_nblocks(); ++i){
        os_H << material.r[1].get_Operator_k()[i] << std::endl;
    }
*/
/*
    H.initialize_fft(*MasterRgrid, material.H.get_Operator_R().get_nrows());
    SettingUp_EigenSystem();
    auto& Uk = Operator<std::complex<double>>::EigenVectors;

    //---------------------------------------------------------------------------------------

    //------------------------setting up TD equations----------------------------------------
    #include "Functional_InitialCondition.hpp"
    #include "Functional_SourceTerm.hpp"
    RK_object.initialize(DensityMatrix, 
                        InitialCondition, SourceTerm);
    RK_object.set_InitialTime(0.);
    RK_object.set_ResolutionTime( ( laser.get_Duration()/laser.get_NumberOfCycles() )/ 1000. );
    PrintResolution = 10;
    kgradient.initialize(*(DensityMatrix.get_Operator(R).get_MeshGrid()));

    coulomb.initialize(material.H.get_Operator_R().get_nrows(), DensityMatrix.get_Operator(R).get_MeshGrid(), material.r);


    print_recap();
    //---------------------------------------------------------------------------------------

    //print DM in R to prove it decays and is zero for large R
    std::stringstream rank;
#ifdef NEGF_MPI
    rank << "DM" << mpi::Communicator::world().rank() << ".txt";
#else
    rank << "DM.txt";
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
*/


template <typename Scalar_T>
void SumWithProduct(Operator<std::complex<double>>& Output, 
                    const Scalar_T& FirstScalar, 
                    const Operator<std::complex<double>>& FirstAddend, 
                    const Scalar_T& SecondScalar, 
                    const Operator<std::complex<double>>& SecondAddend)
{
    PROFILE("SumWithProduct");
    auto& SpaceOfPropagation = Operator<std::complex<double>>::SpaceOfPropagation;

    //std::cout << "space of propagation: "<< (SpaceOfPropagation==k? "k":"R") << std::endl;
    auto& Output_ = Output.get_Operator( SpaceOfPropagation );
    auto& FirstAddend_ = FirstAddend.get_Operator( SpaceOfPropagation );
    auto& SecondAddend_ = SecondAddend.get_Operator( SpaceOfPropagation );

    SumWithProduct(Output_, FirstScalar, FirstAddend_, SecondScalar, SecondAddend_);
}