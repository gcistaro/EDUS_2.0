template<class T>
Simulation::Simulation(const std::string& FileName, const T& arg_meshinit)
{
    PROFILE("Simulation::Initialize");
    //---------------------------getting info from tb----------------------------------------
    material = Material(FileName);
    //---------------------------------------------------------------------------------------
    
    //------------------------------initializing laser---------------------------------------
    laser.set_Intensity(1.e+05, Wcm2);//1.e+16, Wcm2);
    laser.set_InitialTime(0., FemtoSeconds);
    laser.set_Lambda(3000, NanoMeters);
    laser.set_NumberOfCycles(2);
    laser.set_Polarization(Coordinate(1,0,0));
    //---------------------------------------------------------------------------------------

    //--------------------------initializing grids and arrays--------------------------------
    auto MasterRgrid = std::make_shared<MeshGrid>(R, arg_meshinit);
    DensityMatrix.initialize_fft(*MasterRgrid, material.H.get_Operator_R().get_nrows());
    H = DensityMatrix; //to initialize dimensions
    SettingUp_EigenSystem();
    Calculate_Velocity();
    auto& Uk = Operator<std::complex<double>>::EigenVectors;
    auto& UkDagger = Operator<std::complex<double>>::EigenVectors_dagger;
    //---------------------------------------------------------------------------------------

    //------------------------setting up TD equations----------------------------------------
    #include "Functional_InitialCondition.hpp"
    #include "Functional_SourceTerm.hpp"
    RK_object.initialize(DensityMatrix, 
                        InitialCondition, SourceTerm);
    RK_object.set_InitialTime(0.);
    RK_object.set_ResolutionTime(0.01);
    PrintResolution = 100;
    kgradient.initialize();

    DensityMatrix.set_SpaceOfPropagation(SpaceOfPropagation);
    
    print_recap();
    //---------------------------------------------------------------------------------------

    //print DM in R to prove it decays and is zero for large R
    std::ofstream os;
    os.open("DM.txt");
    for(int i=0; i< DensityMatrix.get_Operator_R().get_MeshGrid()->get_mesh().size(); i++){
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
    auto SpaceOfPropagation = Output.get_SpaceOfPropagation();
    auto& Output_ = Output.get_Operator( SpaceOfPropagation );
    auto& FirstAddend_ = FirstAddend.get_Operator( SpaceOfPropagation );
    auto& SecondAddend_ = SecondAddend.get_Operator( SpaceOfPropagation );

    SumWithProduct(Output_, FirstScalar, FirstAddend_, SecondScalar, SecondAddend_);
}