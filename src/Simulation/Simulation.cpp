#include "Simulation/Simulation.hpp"
#include "core/mpi/Communicator.hpp"
#include "core/projectdir.hpp"
#include <cstdlib>
#include <filesystem>

Simulation::Simulation(std::shared_ptr<Simulation_parameters>& ctx__)
{
    PROFILE("Simulation::Initialize");

    output::print("-> initializing simulation");
    ctx_ = ctx__;

    /* set times in a.u. */
    ctx_->cfg().initialtime(Convert(ctx_->cfg().initialtime(), unit(ctx_->cfg().initialtime_units()),
        AuTime));
    ctx_->cfg().finaltime(Convert(ctx_->cfg().finaltime(), unit(ctx_->cfg().finaltime_units()),
        AuTime));
    ctx_->cfg().dt(Convert(ctx_->cfg().dt(), unit(ctx_->cfg().dt_units()),
        AuTime));
    ctx_->cfg().initialtime_units("autime");
    ctx_->cfg().finaltime_units("autime");
    ctx_->cfg().dt_units("autime");

    if ( ctx_->cfg().printresolution_pulse() == 0 ) {
        ctx_->cfg().printresolution_pulse(ctx_->cfg().printresolution());
    }

    output::print("-> initializing material");
    material_ = Material(ctx_->cfg().tb_file());

    output::print("-> initializing grid and arrays");
    auto MasterRgrid = std::make_shared<MeshGrid>(R, ctx_->cfg().grid());
    coulomb_.set_DoCoulomb(ctx_->cfg().coulomb());

    /* getting rytova keldysh with python */
    // ==if (ctx_->cfg().coulomb()) {
    // ==    std::stringstream command;
    // ==    auto grid = ctx_->cfg().grid();
    // ==    auto epsilon = ctx_->cfg().epsilon();
    // ==    auto r0 = ctx_->cfg().r0();
    // ==    command << "python3 " << ProjectDirectory << "/Postproces/RytovaKeldysh.py ";
    // ==    command << grid[0] << " " << grid[1] << " " << grid[2] << " ";
    // ==    command << ctx_->cfg().tb_file() << "_tb.dat";
    // ==    command << " " << r0 << " " << epsilon;
    // ==    output::print("-> creating Rytova Keldish file with python");
    // ==    //output::print(command.str());
    // ==    system(command.str().c_str());
    // ==}
    Operator<std::complex<double>>::SpaceOfPropagation = SpaceOfPropagation_;
    auto& HR = material_.H.get_Operator_R();
    Operator<std::complex<double>>::mpindex.initialize(MasterRgrid->get_Size(), HR.get_nrows() * HR.get_nrows());

    output::print("-> Initializing fft and dft");
    DensityMatrix_.initialize_fft(*MasterRgrid, HR.get_nrows());
    material_.H.dft(DensityMatrix_.get_FT_meshgrid_k().get_mesh(), +1);
    for (auto ix : { 0, 1, 2 }) {
        material_.r[ix].dft(DensityMatrix_.get_FT_meshgrid_k().get_mesh(), +1);
        material_.r[ix].get_Operator(Space::k).make_hermitian();
    }

    output::print("-> initialize fft");
    H_.initialize_fft(DensityMatrix_);

    output::print("-> initializing lasers");
    for (int ilaser = 0; ilaser < int(ctx_->cfg().lasers().size()); ++ilaser) {
        auto currentdata = ctx_->cfg().lasers(ilaser);
        Laser laser;
        
        laser.set_InitialTime(currentdata.t0(), unit(currentdata.t0_units()));
        laser.set_Intensity(currentdata.intensity(), unit(currentdata.intensity_units()));
        auto freq_wavelength = wavelength_or_frequency(ilaser);
        if (freq_wavelength == "frequency") {
            laser.set_Omega(currentdata.frequency(),
                unit(currentdata.frequency_units()));
            ctx_->cfg().dict()["lasers"][ilaser]["wavelength"] = laser.get_Lambda();
        } else {
            laser.set_Lambda(currentdata.wavelength(),
                unit(currentdata.wavelength_units()));
            ctx_->cfg().dict()["lasers"][ilaser]["frequency"] = laser.get_Omega();
        }
        laser.set_NumberOfCycles(currentdata.cycles());
        Coordinate pol(currentdata.polarization()[0], currentdata.polarization()[1],
            currentdata.polarization()[2]);
        pol = pol / pol.norm();
        laser.set_Polarization(pol);
        laser.print_info();
        setoflaser_.push_back(laser);
    }

    output::print("-> solve eigensystem");
    SettingUp_EigenSystem();
    auto& Uk = Operator<std::complex<double>>::EigenVectors;

/* setting up TD equations */
#include "Functional_InitialCondition.hpp"
#include "Functional_SourceTerm.hpp"
    DEsolver_DM_.initialize(DensityMatrix_,
        InitialCondition, SourceTerm,
        solver.at(ctx_->cfg().solver()),
        ctx_->cfg().order());
    DEsolver_DM_.set_InitialTime(Convert(ctx_->cfg().initialtime(),
        unit(ctx_->cfg().initialtime_units()),
        AuTime));
    DEsolver_DM_.set_ResolutionTime(Convert(ctx_->cfg().dt(),
        unit(ctx_->cfg().dt_units()),
        AuTime));

    kgradient_.initialize(*(DensityMatrix_.get_Operator(SpaceOfPropagation_Gradient_).get_MeshGrid()));
    coulomb_.initialize(material_.H.get_Operator_R().get_nrows(),
        DensityMatrix_.get_Operator(R).get_MeshGrid(),
        material_.r);
    DensityMatrix_.go_to_R();
    coulomb_.set_DM0(DensityMatrix_);

    print_recap();

    /* create output directory and cd there */
    std::filesystem::create_directories(std::string(std::filesystem::current_path()) + "/Output");
    std::filesystem::current_path(std::string(std::filesystem::current_path()) + "/Output");

    /* print DM in real space */
    DensityMatrix_.print_Rdecay("DM", material_.rwann_);

    Calculate_Velocity();
    DensityMatrix_.go_to_k();
    assert(DensityMatrix_.get_Operator_k().is_hermitian());

    /* open files for small outputs */
    os_Pop_.open("Population.txt");
    os_Laser_.open("Laser.txt");
    os_VectorPot_.open("Laser_A.txt");
    os_Time_.open("Time.txt");
    os_Velocity_.open("Velocity.txt");

    // ---------------------- create recap file ---------------------------------------------
    ctx_->cfg().dict()["num_bands"] = DensityMatrix_.get_Operator_R().get_nrows();
    ctx_->cfg().dict()["num_kpoints"] = DensityMatrix_.get_Operator_k().get_MeshGrid()->get_TotalSize();

    mdarray<double, 2> bare_k({ DensityMatrix_.get_Operator_k().get_MeshGrid()->get_TotalSize(), 3 });
    for (int ik = 0; ik < DensityMatrix_.get_Operator_k().get_MeshGrid()->get_mesh().size(); ++ik) {
        for (auto& ix : { 0, 1, 2 }) {
            bare_k(ik, ix) = (*(DensityMatrix_.get_Operator_k().get_MeshGrid()))[ik].get(LatticeVectors(k))[ix];
        }
    }
    ctx_->cfg().dict()["kpoints"] = bare_k;
    ctx_->cfg().dict()["A"] = Coordinate::get_Basis(LatticeVectors(R)).get_M();
    ctx_->cfg().dict()["B"] = Coordinate::get_Basis(LatticeVectors(k)).get_M();

#ifdef EDUS_HDF5
    std::string name = "output.h5";
    if (!file_exists(name)) {
        HDF5_tree(name, hdf5_access_t::truncate);
    }
    HDF5_tree fout(name, hdf5_access_t::truncate);
    dump_json_in_h5(ctx_->cfg().dict(), name);
    fout.create_node(nodename::fullH);
    fout.create_node(nodename::DMk);
    fout.create_node(nodename::DMk_bloch);
    fout.create_node(nodename::SelfEnergy);
    fout.create_node(nodename::time_au);
    mpi::Communicator::world().barrier();
#endif
    //---------------------------------------------------------------------------------------
}

bool Simulation::PrintObservables(const double& time__, const bool& use_sparse)
{
    /* check if we are within (any) pulse */
    int printresolution;
    if( !use_sparse ) {
        printresolution = ctx_->cfg().printresolution_pulse();
    } else {
        for (int ilaser = 0; ilaser < setoflaser_.size(); ++ilaser) {
            if (time__ > setoflaser_[ilaser].get_InitialTime() + 1.e-07 && time__ < setoflaser_[ilaser].get_FinalTime() + 1.e-07) {
                printresolution = ctx_->cfg().printresolution_pulse();
                break;
            } else {
                printresolution = ctx_->cfg().printresolution();
            }
        }
    }
    return (int(round(time__ / DEsolver_DM_.get_ResolutionTime())) % printresolution == 0);
}

void Simulation::SettingUp_EigenSystem()
{
    PROFILE("Simulation::SetEigensystem");

    //------------------------go to H(k)--------------------------------------------
    auto& MasterkGrid = DensityMatrix_.get_Operator_k().get_MeshGrid()->get_mesh();
    // material.H.dft(MasterkGrid, +1);
    //------------------------------------------------------------------------------

    //--------------------solve eigen problem---------------------------------------
    auto& Uk = Operator<std::complex<double>>::EigenVectors;
    auto& UkDagger = Operator<std::complex<double>>::EigenVectors_dagger;

    material_.H.get_Operator_k().diagonalize(Band_energies_, Uk);

    /* open the gap */
    


    //------------------------------------------------------------------------------

    //-------------------get U dagger-----------------------------------------------
    UkDagger.initialize(k, Uk.get_nblocks(), Uk.get_ncols(), Uk.get_nrows());
    for (int ik = 0; ik < UkDagger.get_nblocks(); ++ik) {
        for (int ir = 0; ir < UkDagger.get_nrows(); ++ir) {
            for (int ic = 0; ic < UkDagger.get_ncols(); ++ic) {
                UkDagger[ik](ir, ic) = conj(Uk[ik](ic, ir));
            }
        }
    }
    //------------------------------------------------------------------------------
};

void Simulation::Calculate_TDHamiltonian(const double& time__, const bool& erase_H__)
{
#ifdef EDUS_TIMERS
    PROFILE("SourceTerm::Calculate_TDHamiltonian");
#endif
    H_.go_to_k();
    //--------------------get aliases for nested variables--------------------------------
    auto& H = H_.get_Operator(SpaceOfPropagation_);
    auto& H0 = material_.H.get_Operator(SpaceOfPropagation_);
    auto& x = material_.r[0].get_Operator(SpaceOfPropagation_);
    auto& y = material_.r[1].get_Operator(SpaceOfPropagation_);
    auto& z = material_.r[2].get_Operator(SpaceOfPropagation_);
    auto las = setoflaser_(time__).get("Cartesian");

    // auto& ci = MeshGrid::ConvolutionIndex[{H0_.get_MeshGrid()->get_id(),
    //                                           H_.get_MeshGrid()->get_id(),
    //                                           Operator<std::complex<double>>::MeshGrid_Null->get_id()}];
    //-----------------------------------------------------------------------------------

    //--------------------------do initializations---------------------------------------
    if (erase_H__) {
        H.fill(0.);
    }

    // if(ci.get_Size(0) == 0 && SpaceOfPropagation == R) {
    //     MeshGrid::Calculate_ConvolutionIndex( *(H0_.get_MeshGrid()),
    //                                               *(H_.get_MeshGrid()),
    //                                               *(Operator<std::complex<double>>::MeshGrid_Null));
    // }
    //---------------------------------------------------------------------------------

    //------------------------H(R) = H0(R) + E.r(R)-------------------------------------
    // == /* H_ = H0_ + las_x \cdot x */
    // == SumWithProduct(H_, 1., H0_, las[0], x_);
    // == /* H_ = H_ + las_y \cdot y */
    // == SumWithProduct(H_, 1., H_, las[1], y_);
    // == /* H_ = H_ + las_z \cdot z */
    // == SumWithProduct(H_, 1., H0_, las[1], y_);

#pragma omp parallel for schedule(static) collapse(3)
    for (int iblock = 0; iblock < H0.get_nblocks(); ++iblock) {
        for (int irow = 0; irow < H0.get_nrows(); ++irow) {
            for (int icol = 0; icol < H0.get_ncols(); ++icol) {
                // auto Hblock = ( SpaceOfPropagation == k ? iblock : ci(iblock, 0) );
                // H_(Hblock, irow, icol) = H0_(iblock, irow, icol)
                H(iblock, irow, icol) += H0(iblock, irow, icol)
                    + las[0] * x(iblock, irow, icol)
                    + las[1] * y(iblock, irow, icol)
                    + las[2] * z(iblock, irow, icol);
            }
        }
    }
    //---------------------------------------------------------------------------------
}

void Simulation::Propagate()
{
    PROFILE("Simulation::Propagate");
    int iFinalTime = int((ctx_->cfg().finaltime() - ctx_->cfg().initialtime()) / DEsolver_DM_.get_ResolutionTime()) + 2;

    /* print info on output */
    output::title("PROPAGATION");
    output::print("Initial time:       *", ctx_->cfg().initialtime(), " a.u.", Convert(ctx_->cfg().initialtime(), AuTime, FemtoSeconds), " fs");
    output::print("Final time:         *", ctx_->cfg().finaltime(), " a.u.", Convert(ctx_->cfg().finaltime(), AuTime, FemtoSeconds), " fs");
    output::print("Number of steps:    *", iFinalTime);

    /* do steps */
    for (int it = 0; it < iFinalTime; ++it) {
        if (it % 100 == 0) {
            output::print("it:                 *", it, " / ", iFinalTime, 100 * double(it) / iFinalTime, " %");
        }
        do_onestep();
    }

    output::print(std::string(output::linesize - 5, '*'));
}

void Simulation::do_onestep()
{
    auto CurrentTime = DEsolver_DM_.get_CurrentTime();
    //------------------------Print population-------------------------------------

    if (PrintObservables(CurrentTime, true)) {
#ifdef EDUS_HDF5
        std::string name = "output.h5";
        std::stringstream node;
        node << get_it_sparse(CurrentTime);

        HDF5_tree fout(name, hdf5_access_t::read_write);

        if(ctx_->cfg().dict()["toprint"]["DMk_wannier"] == "true") {
            DensityMatrix_.go_to_k(false);
            DensityMatrix_.get_Operator_k().write_h5(name, nodename::DMk, node.str());
        }
        if(ctx_->cfg().dict()["toprint"]["DMk_bloch"] == "true") {
            DensityMatrix_.go_to_k(false);
            DensityMatrix_.go_to_bloch();
            DensityMatrix_.get_Operator_k().write_h5(name, nodename::DMk_bloch, node.str());
            DensityMatrix_.go_to_wannier();
        } 
        if(ctx_->cfg().dict()["toprint"]["SelfEnergy"] == "true") {
            DensityMatrix_.go_to_R(true);
            H_.go_to_R(false);
            coulomb_.EffectiveHamiltonian(H_, DensityMatrix_, true);
            H_.go_to_k(true);
            H_.get_Operator_k().write_h5(name, nodename::SelfEnergy, node.str());
            DensityMatrix_.go_to_k(false);
        }
        if(ctx_->cfg().dict()["toprint"]["fullH"] == "true") {
            Calculate_TDHamiltonian(CurrentTime, true);
            H_.get_Operator_k().write_h5(name, nodename::fullH, node.str());
        }
        fout[nodename::time_au].write(node.str(), CurrentTime);
#endif
    }

    if (PrintObservables(CurrentTime, false)) {
        // print time
        os_Time_ << CurrentTime << std::endl;
        // print laser
        os_Laser_ << setoflaser_(DEsolver_DM_.get_CurrentTime()).get("Cartesian");
        os_VectorPot_ << setoflaser_.VectorPotential(DEsolver_DM_.get_CurrentTime()).get("Cartesian");
        Print_Population();
        Print_Velocity();
    }
    //------------------------------------------------------------------------------
    DEsolver_DM_.Propagate();
}

void Simulation::Print_Population()
{
    static Operator<std::complex<double>> aux_DM;
    aux_DM = DensityMatrix_;

    aux_DM.go_to_bloch();

    auto Population = TraceK(aux_DM.get_Operator(Space::k));
#ifdef EDUS_MPI
    if (kpool_comm.rank() == 0)
#endif
    {
        for (int ibnd = 0; ibnd < int(Population.size()); ibnd++) {
            Population[ibnd] /= double(aux_DM.get_Operator_k().get_MeshGrid()->get_TotalSize());
        }
        for (int ibnd = 0; ibnd < int(Population.size()); ibnd++) {
            if (ibnd < ctx_->cfg().filledbands())
                Population[ibnd] = 1. - Population[ibnd];
            os_Pop_ << std::setw(20) << std::setprecision(6) << Population[ibnd].real();
            os_Pop_ << " ";
        }
        os_Pop_ << std::endl;
    }
}


double Simulation::jacobian(const Matrix<double>& A__) const
{
    /* check if we are considering a slab or a 3D material, from the k grid */
    double j = 0.0;
    auto dim = 0;
    for (auto& ix : {0,1,2}) {
        dim += ( ctx_->cfg().grid()[ix] > 1 ? 1 : 0 ); 
    }

    if ( dim == 3 ) {
        j = std::abs(A__.determinant());
    }
    else if ( dim == 2 ) {
        /* check which lattice vectors is perpendicular to the slab */
        int perpendicular = -1;
        for (auto& ix : {0,1,2}) {
            if( ctx_->cfg().grid()[ix] <= 1 ) {
                perpendicular = ix; 
                break;
            }
        }
        if(perpendicular != 2) {
            std::stringstream ss;
            ss << "perpendicular = " << perpendicular << " !=2 -> implement me!\n";
            throw std::runtime_error(ss.str());
        }
        j = std::abs(A__(0,0)*A__(1,1) - A__(1,0)*A__(0,1));
    }
    else if ( dim == 1 ) {
        throw std::runtime_error("Implement me!\n");
    }
    return j;
}

void Simulation::Calculate_Velocity()
{
    Calculate_TDHamiltonian(-6000, true);
    H_.go_to_R();
    std::vector<Coordinate> direction(3);
    direction[0].initialize(1, 0, 0);
    direction[1].initialize(0, 1, 0);
    direction[2].initialize(0, 0, 1);

    for (int ix : { 0, 1, 2 }) {
        Velocity_[ix].initialize_fft(DensityMatrix_);
        Velocity_[ix].lock_space(k);
        Velocity_[ix].get_Operator_k().fill(0.);

        /* V = -i*[r,H0] */
        commutator(Velocity_[ix].get_Operator_k(), -im, material_.r[ix].get_Operator_k(), H_.get_Operator_k());

        if (SpaceOfPropagation_Gradient_ == R ) {
            Velocity_[ix].go_to_R();
        }
        /* V += R*H0 */
        kgradient_.Calculate(1., Velocity_[ix].get_Operator(SpaceOfPropagation_Gradient_), 
                                H_.get_Operator(SpaceOfPropagation_Gradient_), direction[ix], false);
        if (SpaceOfPropagation_Gradient_ == R ) {
            Velocity_[ix].go_to_k();
        }

        /* apply right constants for the normalization of wavefunctions */
        auto det = jacobian(Coordinate::get_Basis(LatticeVectors(Space::R)).get_M());
        #pragma omp parallel for 
        for (int iblock = 0; iblock < Velocity_[ix].get_Operator(Space::k).get_nblocks(); ++iblock) {
            for (int irow = 0; irow < Velocity_[ix].get_Operator(Space::k).get_nrows(); ++irow) {
                for (int icol = 0; icol < Velocity_[ix].get_Operator(Space::k).get_ncols(); ++icol) {
                    Velocity_[ix].get_Operator(Space::k)(iblock, irow, icol) /= det;
                }
            }
        }   
    }
}

void Simulation::Print_Velocity()
{
    std::array<std::complex<double>, 3> v = { 0., 0., 0. };

    auto& DMK = DensityMatrix_.get_Operator_k();
    static BlockMatrix<std::complex<double>> temp(k, DMK.get_nblocks(), DMK.get_nrows(), DMK.get_ncols());

    for (auto ix : { 0, 1, 2 }) {
        temp.fill(0.);
        multiply(temp, 1. + im * 0., Velocity_[ix].get_Operator_k(), DMK);
        auto ToSum = TraceK(temp);
        for (int ibnd = 0; ibnd < int(ToSum.size()); ++ibnd) {
            v[ix] += ToSum[ibnd];
        }
        v[ix] /= DMK.get_MeshGrid()->get_TotalSize();
    }
#ifdef EDUS_MPI
    if (kpool_comm.rank() == 0)
#endif
    {
        os_Velocity_ << std::setw(20) << std::setprecision(8) << v[0].real();
        os_Velocity_ << std::setw(20) << std::setprecision(8) << v[0].imag();
        os_Velocity_ << std::setw(20) << std::setprecision(8) << v[1].real();
        os_Velocity_ << std::setw(20) << std::setprecision(8) << v[1].imag();
        os_Velocity_ << std::setw(20) << std::setprecision(8) << v[2].real();
        os_Velocity_ << std::setw(20) << std::setprecision(8) << v[2].imag();
        os_Velocity_ << std::endl;
    }
}

void Simulation::print_recap()
{
    std::stringstream title;
    title << std::string(5, ' ') << "INPUT RECAP" << std::string(5, ' ');
    int num_stars = (output::linesize - title.str().length()) / 2 - 2;

    output::title("INPUT RECAP");
    //==    output::print("input file:         *", std::string(8, ' '), JsonFile_ );
    output::print("tb_model                 *", std::string(8, ' '), ctx_->cfg().tb_file());
    output::print("grid                     *", std::string(8, ' '), "[",
        DensityMatrix_.get_Operator_R().get_MeshGrid()->get_Size()[0], ", ",
        DensityMatrix_.get_Operator_R().get_MeshGrid()->get_Size()[1], ", ",
        DensityMatrix_.get_Operator_R().get_MeshGrid()->get_Size()[2], "]");
    output::print("Resolution time          *", DEsolver_DM_.get_ResolutionTime(), " a.u.",
        Convert(DEsolver_DM_.get_ResolutionTime(), AuTime, FemtoSeconds), " fs");
    output::print("PrintResolution          *", ctx_->cfg().printresolution());
    output::print("PrintResolution(pulse):  *", ctx_->cfg().printresolution_pulse());
    output::print("Coulomb                  *", std::string(8, ' '), (coulomb_.get_DoCoulomb() ? "True" : "False"));
    output::print("epsilon                  *", ctx_->cfg().epsilon());
    output::print("r0                       *", Convert(ctx_->cfg().r0(), Angstrom, AuLength), " a.u.",
        ctx_->cfg().r0(), " angstrom");
    output::print("toprint-> DMk_wannier    *        ", std::string(ctx_->cfg().dict()["toprint"]["DMk_wannier"]));
    output::print("toprint-> DMk_bloch      *        ", std::string(ctx_->cfg().dict()["toprint"]["DMk_bloch"]));
    output::print("toprint-> fullH          *        ", std::string(ctx_->cfg().dict()["toprint"]["fullH"]));
    output::print("toprint-> SelfEnergy     *        ", std::string(ctx_->cfg().dict()["toprint"]["SelfEnergy"]));
    output::stars();

    output::stars();
    output::title("WANNIER");

    output::print("#R points :         *", material_.H.get_Operator_R().get_nblocks());
    auto& A = Coordinate::get_Basis(LatticeVectors(R)).get_M();
    output::print("          :         *", A(0, 0), A(0, 1), A(0, 2));
    output::print("    A     :         *", A(1, 0), A(1, 1), A(1, 2));
    output::print("          :         *", A(2, 0), A(2, 1), A(2, 2));
    auto& B = Coordinate::get_Basis(LatticeVectors(k)).get_M();
    output::print("          :         *", B(0, 0), B(0, 1), B(0, 2));
    output::print("    B     :         *", B(1, 0), B(1, 1), B(1, 2));
    output::print("          :         *", B(2, 0), B(2, 1), B(2, 2));
    output::print("");
    output::print(" jacobian(A)        *", jacobian(A));
    output::print("");
    output::stars();

    for (int ilaser = 0; ilaser < int(setoflaser_.size()); ++ilaser) {
        setoflaser_[ilaser].print_info();
    }
}

int Simulation::get_it(const double& time__) const
{
    return int(time__ / DEsolver_DM_.get_ResolutionTime());
}

int Simulation::get_it_sparse(const double& time__) const
{
    static int counter = -1;
    counter++;
    return counter;
    // return int(round(time__/DEsolver_DM_.get_ResolutionTime()/ctx_->cfg().printresolution()));
}

void Simulation::print_grids()
{
    auto& MasterRgrid = DensityMatrix_.get_Operator_R().get_MeshGrid();
    auto& kgrid = DensityMatrix_.get_Operator_k().get_MeshGrid();
    auto Rfft = std::make_shared<MeshGrid>(fftPair(*kgrid)); // Rspace as fourier space of kgrid
    std::ofstream os;
    os.open("MasterR.txt");
    for (auto& R : MasterRgrid->get_mesh()) {
        os << R.get("Cartesian");
    }
    os.close();
    os.open("k.txt");
    for (auto& k : kgrid->get_mesh()) {
        os << k.get("Cartesian");
    }
    os.close();
    auto MeshRfft = fftPair(*kgrid);
    os.open("R.txt");
    for (auto& R : Rfft->get_mesh()) {
        os << R.get("Cartesian");
    }
    os.close();
}

std::string Simulation::wavelength_or_frequency(const int& idx__)
{
    bool is_frequency = (std::abs(ctx_->cfg().lasers(idx__).frequency()) > 1.e-07);
    if (is_frequency)
        return "frequency";
    bool is_wavelength = (ctx_->cfg().lasers(idx__).wavelength());
    if (!is_wavelength) {
        throw std::runtime_error("You must specify (nonzero) frequency *xor* wavelength!");
    }
    return "wavelength";
}
