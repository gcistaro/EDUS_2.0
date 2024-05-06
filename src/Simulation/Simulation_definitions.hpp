template<class T>
Simulation::Simulation(const std::string& FileName, const T& arg_meshinit)
{
    PROFILE_START("Simulation::Initialize");
    material = Material(FileName);
    
    laser.set_Intensity(1.e+05, Wcm2);//1.e+16, Wcm2);
    laser.set_InitialTime(0., FemtoSeconds);
    laser.set_Lambda(3000, NanoMeters);
    laser.set_NumberOfCycles(2);
    laser.set_Polarization(Coordinate<k>(1,0,0));

    //auto MasterRgrid = std::make_shared<MeshGrid<R>>(Radius);//spherical grid for exponentially decaying DM
    //auto MasterRgrid = std::make_shared<MeshGrid<R>>(std::array<int,3>{30,30,1});//spherical grid for exponentially decaying DM
    auto MasterRgrid = std::make_shared<MeshGrid<R>>(arg_meshinit);

    DensityMatrix.initialize_fft(*MasterRgrid, material.H.get_Operator_R().get_nrows());

    H = DensityMatrix; //to initialize dimensions
    //Calculate_Position();

    SettingUp_EigenSystem();
    Calculate_Velocity();
    auto& Uk = Operator<std::complex<double>>::EigenVectors;
    auto& UkDagger = Operator<std::complex<double>>::EigenVectors_dagger;
   
    #include "Functional_InitialCondition.hpp"
    #include "Functional_SourceTerm.hpp"
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
    RK_object.set_ResolutionTime(0.001);
    print_recap();
    PROFILE_STOP("Simulation::Initialize");
    //this->Propagate();
}
