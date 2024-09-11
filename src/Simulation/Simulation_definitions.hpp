template<class T>
Simulation::Simulation(const std::string& FileName, const T& arg_meshinit)
{
    PROFILE("Simulation::Initialize");
    //---------------------------getting info from tb----------------------------------------
    material = Material(FileName);
    //---------------------------------------------------------------------------------------
    
    //------------------------------initializing laser---------------------------------------
    laser.set_Intensity(1.e+11, Wcm2);//1.e+16, Wcm2);
    laser.set_InitialTime(0., FemtoSeconds);
    laser.set_Intensity(1.e+05, Wcm2);
    laser.set_Lambda(800, NanoMeters);
    laser.set_NumberOfCycles(10);    
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
    H.initialize_dims(DensityMatrix);
    SettingUp_EigenSystem();
    //Calculate_Velocity();
    auto& Uk = Operator<std::complex<double>>::EigenVectors;

    //---------------------------------------------------------------------------------------

    //------------------------setting up TD equations----------------------------------------
    #include "Functional_InitialCondition.hpp"
    #include "Functional_SourceTerm.hpp"
    RK_object.initialize(DensityMatrix, 
                        InitialCondition, SourceTerm);
    RK_object.set_InitialTime(0.);
    RK_object.set_ResolutionTime(0.01);
    PrintResolution = 1;
    kgradient.initialize(*(DensityMatrix.get_Operator(R).get_MeshGrid()));

    

    print_recap();
    //---------------------------------------------------------------------------------------

    //print DM in R to prove it decays and is zero for large R
    std::ofstream os;
    os.open("DM.txt");
    for(int i=0; i< DensityMatrix.get_Operator_R().get_nblocks(); ++i){
        os << DensityMatrix.get_Operator_R().get_MeshGrid()->get_mesh()[i].norm();
        os << " " << std::abs(max(DensityMatrix.get_Operator_R()[i])) << std::endl;
    }
    os.close();
    //---------------------------------------------------------------------------------------

    //---------------------open files---------------------------------------------------------
    os_Pop.open("Population.txt");
    os_Laser.open("Laser.txt");
    os_Time.open("Time.txt");
    //---------------------------------------------------------------------------------------
}



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