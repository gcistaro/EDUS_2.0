#include "Simulation/Simulation.hpp"
#include "core/mpi/Communicator.hpp"
#include "core/projectdir.hpp"
#include <cstdlib>
#include <filesystem>

/// @brief Initialization of all the variables needed for the simulation, as well as allocation of memory
/// and definition of the differential equations
/// @param ctx__ Simulation parameters already set up and ready to initialize everything in the simulation class
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

    /* set r0 in a.u. */
    std::vector<double> r0_au(ctx_->cfg().r0().size());
    for ( int i=0; i<ctx_->cfg().r0().size(); ++i ) {
        r0_au[i] = Convert(ctx_->cfg().r0()[i], unit(ctx_->cfg().r0_units()), AuLength);
    }
    ctx_->cfg().r0(r0_au);

    if ( ctx_->cfg().printresolution_pulse() == 0 ) {
        ctx_->cfg().printresolution_pulse(ctx_->cfg().printresolution());
    }

    ctx_->cfg().opengap(Convert(ctx_->cfg().opengap(), unit(ctx_->cfg().opengap_units()),
        AuEnergy));
    ctx_->cfg().opengap_units("auenergy");

    SpaceOfPropagation_Gradient_ = (ctx_->cfg().gradient_space() == "R" ? Space::R : Space::k);

    output::print("-> initializing material");
    material_ = Material(ctx_->cfg().tb_file());

    output::print("-> initializing grid and arrays");
    MeshGrid::MasterRgrid = MeshGrid(R, ctx_->cfg().grid());
    MeshGrid::MasterRgrid_GammaCentered = get_GammaCentered_grid(MeshGrid::MasterRgrid);
    coulomb_.set_DoCoulomb(ctx_->cfg().coulomb());
    coulomb_.set_epsilon(ctx_->cfg().epsilon());
    coulomb_.set_r0(ctx_->cfg().r0());

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
    Operator<std::complex<double>>::mpindex.initialize(MeshGrid::MasterRgrid.get_Size(), HR.get_nrows() * HR.get_nrows());

    output::print("-> Initializing fft and dft");
    DensityMatrix_.initialize_fft(MeshGrid::MasterRgrid, HR.get_nrows());

#ifdef __DEBUG_MODE
    for(int ik=0; ik<DensityMatrix_.get_Operator(Space::k).get_MeshGrid()->get_TotalSize(); ++ik){
        for(int iR=0; iR<DensityMatrix_.get_Operator(Space::R).get_MeshGrid()->get_TotalSize(); ++iR) {
            auto& kp = (*(DensityMatrix_.get_Operator(Space::k).get_MeshGrid()))[ik];
            auto& Rp = (*(DensityMatrix_.get_Operator(Space::R).get_MeshGrid()))[iR];
            auto& kcrys = kp.get(LatticeVectors(k));
            auto& Rcrys = Rp.get(LatticeVectors(R));
            if( std::abs(kp.dot(Rp) - 2.*pi*(kcrys[0]*Rcrys[0]+kcrys[1]*Rcrys[1]+kcrys[2]*Rcrys[2])) > 1.e-14) {
                std::cout << kp.get("Cartesian") << Rp.get("Cartesian");
                std::cout << kcrys << Rcrys;
                std::cout << ik << " " << iR <<" "  << kp.dot(Rp) << " " << 2.*pi*(kcrys[0]*Rcrys[0]+kcrys[1]*Rcrys[1]+kcrys[2]*Rcrys[2]) <<   std::endl;
                std::cout << std::endl;
            }
        }
    }
#endif

    material_.H.dft(DensityMatrix_.get_FT_meshgrid_k().get_mesh(), +1);
    for (auto ix : { 0, 1, 2 }) {
        material_.r[ix].dft(DensityMatrix_.get_FT_meshgrid_k().get_mesh(), +1);
        material_.r[ix].get_Operator(Space::k).make_hermitian();
    }

    H_.initialize_fft(DensityMatrix_);
    H0_.initialize_fft(DensityMatrix_);
    r_[0].initialize_fft(DensityMatrix_);
    r_[1].initialize_fft(DensityMatrix_);
    r_[2].initialize_fft(DensityMatrix_);

    auto& materialH0k = material_.H.get_Operator(Space::k);
    auto& materialr0k  = material_.r[0].get_Operator(Space::k);
    auto& materialr1k  = material_.r[1].get_Operator(Space::k);
    auto& materialr2k  = material_.r[2].get_Operator(Space::k);
    std::copy(materialH0k.begin(), materialH0k.end(), H0_.get_Operator(Space::k).begin());
    std::copy(materialr0k.begin(), materialr0k.end(), r_[0].get_Operator(Space::k).begin());
    std::copy(materialr1k.begin(), materialr1k.end(), r_[1].get_Operator(Space::k).begin());
    std::copy(materialr2k.begin(), materialr2k.end(), r_[2].get_Operator(Space::k).begin());

    H0_.lock_space(Space::k);           H0_.go_to_R();
    r_[0].lock_space(Space::k);         r_[0].go_to_R();
    r_[1].lock_space(Space::k);         r_[1].go_to_R();
    r_[2].lock_space(Space::k);         r_[2].go_to_R();

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
        laser.set_Phase(currentdata.phase());
        Coordinate pol(currentdata.polarization()[0], currentdata.polarization()[1],
            currentdata.polarization()[2]);
        pol = pol / pol.norm();
        laser.set_Polarization(pol);
        setoflaser_.push_back(laser);
    }

    output::print("-> solve eigensystem");
    SettingUp_EigenSystem();
    if( ctx_->cfg().opengap() ) OpenGap();
    auto& Uk = Operator<std::complex<double>>::EigenVectors;
    if( ctx_->cfg().kpath().size() > 1 ) {
        output::print("-> Printing band structure");
        print_bandstructure(ctx_->cfg().kpath(), material_.H);
    }


    if( ctx_->cfg().kpath().size() > 1 ) {
        output::print("-> Printing band structure");
        print_bandstructure(ctx_->cfg().kpath(), material_.H);
    }


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

/// @brief Defines whether or not at the current time step time__ we print the txt files and the matrices in hdf5
/// @param time__ Current time in the simulation
/// @param use_sparse If true, we never use the sparse printing, mainly needed for the printing of big matrices
/// when laser is off
/// @return Boolean defining if we want to print the observables
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

/// @brief Driver for the diagonalization of the Hamiltonian.
/// Here we already know the k grid on which we want to propagate our system, defined in the density matrix.
/// With that and H0(k) already set up in the material class (material_.H), we are able to diagonalize
/// and getting eigenvalues (Band_energies_) and eigenvectors (Uk) that will be used later.
void Simulation::SettingUp_EigenSystem()
{
    PROFILE("Simulation::SetEigensystem");

    auto& MasterkGrid = DensityMatrix_.get_Operator_k().get_MeshGrid()->get_mesh();

    /* solve eigen problem */
    auto& Uk = Operator<std::complex<double>>::EigenVectors;
    auto& UkDagger = Operator<std::complex<double>>::EigenVectors_dagger;

    material_.H.get_Operator_k().diagonalize(Band_energies_, Uk);

    /* get U dagger */
    UkDagger.initialize(k, Uk.get_nblocks(), Uk.get_ncols(), Uk.get_nrows());
    for (int ik = 0; ik < UkDagger.get_nblocks(); ++ik) {
        for (int ir = 0; ir < UkDagger.get_nrows(); ++ir) {
            for (int ic = 0; ic < UkDagger.get_ncols(); ++ic) {
                UkDagger[ik](ir, ic) = conj(Uk[ik](ic, ir));
            }
        }
    }

    /* TODO: open the gap */
};

/// @brief Function to calculate the time dependent Hamiltonian (one-body) as a sum of H0 and the
/// interaction with the laser:
/// @f[
/// H_{\text{1B}}(\textbf{k}) = H_0(\textbf{k})+ \boldsymbol{\varepsilon}(t)\cdot \Xi(\textbf{k})
///@f]
/// @param time__ Time on which we want to calculate the hamiltonian, used for the laser
/// @param erase_H__ If true, we delete H_ before computing it. Needs to be false during
/// time propagation to sum up contributions
void Simulation::Calculate_TDHamiltonian(const double& time__, const bool& erase_H__)
{
#ifdef EDUS_TIMERS
    PROFILE("SourceTerm::Calculate_TDHamiltonian");
#endif
    //--------------------get aliases for nested variables--------------------------------
    auto& H = H_.get_Operator(SpaceOfCalculateTDHamiltonian_);
    auto& H0 = H0_.get_Operator(SpaceOfCalculateTDHamiltonian_);
    auto& x = r_[0].get_Operator(SpaceOfCalculateTDHamiltonian_);
    auto& y = r_[1].get_Operator(SpaceOfCalculateTDHamiltonian_);
    auto& z = r_[2].get_Operator(SpaceOfCalculateTDHamiltonian_);

    auto las  = setoflaser_(time__).get("Cartesian");

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
    // == SumWithProduct(H, 1., H0, las[0], x);
    // == /* H_ = H_ + las_y \cdot y */
    // == SumWithProduct(H, 1., H, las[1], y);
    // == /* H_ = H_ + las_z \cdot z */
    // == SumWithProduct(H, 1., H, las[2], z);

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
    H_.lock_space(SpaceOfCalculateTDHamiltonian_);
}

/// @brief Driver for the full time propagation. It prints every 100 steps in the standard output a message.
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
    double start_time = omp_get_wtime();

    for (int it = 0; it < iFinalTime; ++it) {
        if (it % 100 == 0) {
            output::print("it:                 *", it, " / ", iFinalTime, 100 * double(it) / iFinalTime, " %",
                          "               time: ", omp_get_wtime() - start_time, " sec");
        }
        do_onestep();
    }

    output::print(std::string(output::linesize - 5, '*'));
}

/// @brief Driver for single time step propagation.
/// - checks if something needs to be printed using PrintObservables
/// - prints .txt file and/or h5
/// - Triggers the one step propagation in DEsolver
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
        Print_Velocity(DensityMatrix_);
    }
    //------------------------------------------------------------------------------
    DEsolver_DM_.Propagate();
}

/// @brief Prints in Population.txt the Population for each band:
/// @f[
/// P_n(t) = \frac{1}{N}\sum_{\textbf{k}} \rho_{nn}(\textbf{k})
/// @f]
/// where @f$ \rho_{nn}(\textbf{k}) @f$ is the density matrix in the bloch gauge.
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

/// @brief Calculation of the jacobian of the real lattice vectors
/// @param A__ Matrix containing the real lattice vectors spanned over the columns
/// @return The jacobian of the matrix. Notice that A__ is always 3x3 while we calculate the jacobian
/// using the submatrix when we have a slab material. In that case only the submatrix 2x2 is used.
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

/// @brief Calculation of the matrix elements of the velocity operator on the basis in k.
/// As the density matrix is propagated in wannier gauge, this is the gauge where we calculate it.
/// The equation implemented here is:
/// @f[ \textbf{v}_{nm}(\textbf{k}) = -i[\textbf{\xi}(\textbf{k}), H_0(\textbf{k})] + \nabla_{\textbf{k}} H_0(\textbf{k}) @f]
/// or, in R space:
/// @f[ \textbf{v}_{nm}(\textbf{R}) = -i[\textbf{r}, H_0(\textbf{k})] + i \textbf{R} H_0(\textbf{k}) @f]
/// It must be noticed that because we want to represent the infinite system, when we sum over k an additional constant
/// need to be added, to take into account the elementary volume in the integral:
/// @f[
/// \delta k = \frac{\Omega_{\text{BZ}}}{N_k} = \frac{(2\pi)^d}{N_k \Omega_{\text{UC}}}
/// @f]
/// where @f$\Omega_{\text{BZ}}@f$ and @f$\Omega_{\text{UC}}@f$ are, respectively, the volume(/area) of the Brillouin zone and of
/// the unit cell, while d is the dimension of the system. In this function we multiply by everything but @f$\frac{1}{N_k}@f$
void Simulation::Calculate_Velocity()
{
#ifdef __DEBUG
    H_.print_Rdecay("H0__", material_.rwann_);
    for (int ix : { 0, 1, 2 }) {
        std::stringstream name;
        name << "r0__" << ix ;
        material_.r[ix].print_Rdecay(name.str(), material_.rwann_);
    }
#endif
    std::vector<Coordinate> direction(3);
    direction[0].initialize(1, 0, 0);
    direction[1].initialize(0, 1, 0);
    direction[2].initialize(0, 0, 1);

    for (int ix : { 0, 1, 2 }) {
        Velocity_[ix].initialize_fft(DensityMatrix_);
        Velocity_[ix].lock_space(k);
        Velocity_[ix].get_Operator_k().fill(0.);

        /* V = -i*[r,H0] */
        commutator(Velocity_[ix].get_Operator_k(), -im, r_[ix].get_Operator_k(), H0_.get_Operator_k());

        if (SpaceOfPropagation_Gradient_ == R ) {
            Velocity_[ix].go_to_R();
        }
#ifdef __DEBUG
        std::stringstream vel;
        vel << "velocity__" << ix;
        Velocity_[ix].print_Rdecay(vel.str(), material_.rwann_);
#endif
        /* V += i*R*H0 */
        kgradient_.Calculate(1., Velocity_[ix].get_Operator(SpaceOfPropagation_Gradient_),
                                H0_.get_Operator(SpaceOfPropagation_Gradient_), direction[ix], false);
        if (SpaceOfPropagation_Gradient_ == R ) {
            Velocity_[ix].go_to_k();
        }
#ifdef __DEBUG
        vel.str("");
        vel << "finalvelocity__" << ix;
        Velocity_[ix].print_Rdecay(vel.str(), material_.rwann_);
#endif
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

/// @brief Calculates the mean value of the velocity operator and prints it in Output/Velocity.txt.
/// The equation implemented is:
/// \textbf{v}(t) = \text{Tr}(\rho \textbf{v})
/// where v are the matrix elements of the velocity operator in the basis.
/// Everything is conveniently done in k space.
void Simulation::Print_Velocity(Operator<std::complex<double>>& aux_DM)
{
    std::array<std::complex<double>, 3> v = { 0., 0., 0. };

    auto& DMR = aux_DM.get_Operator(Space::R);
    static mdarray<std::complex<double>,1> Peierls_phase({DMR.get_nblocks()});
// ==     if ( ctx_->cfg().peierls() ) {
// ==         aux_DM.go_to_R();
// ==         auto& Rgrid = *(DMR.get_MeshGrid());
// ==         auto lasA = setoflaser_.VectorPotential(DEsolver_DM_.get_CurrentTime());
// == #pragma omp parallel for schedule(static)
// ==         for (int iR_loc = 0; iR_loc < DMR.get_nblocks(); ++iR_loc) {
// ==             int iR_glob = Rgrid.mpindex.loc1D_to_glob1D(iR_loc);
// ==             Peierls_phase(iR_loc) = std::exp(-im*lasA.dot(Rgrid[iR_glob]));
// ==         }
// == 
// == #pragma omp parallel for schedule(static) collapse(3)
// ==         for (int iblock = 0; iblock < DMR.get_nblocks(); ++iblock) {
// ==             for (int irow = 0; irow < DMR.get_nrows(); ++irow) {
// ==                 for (int icol = 0; icol < DMR.get_ncols(); ++icol) {
// ==                     DMR(iblock, irow, icol) *= Peierls_phase(iblock);
// ==                 }
// ==             }
// ==         }
// ==         aux_DM.go_to_k();
// ==     }
    auto& DMK =  aux_DM.get_Operator(Space::k);

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
//==     if ( ctx_->cfg().peierls() ) {
//==         aux_DM.go_to_R();
//==         auto& Rgrid = *(DMR.get_MeshGrid());
//==         auto lasA = setoflaser_.VectorPotential(DEsolver_DM_.get_CurrentTime());
//== #pragma omp parallel for schedule(static)
//==         for (int iR_loc = 0; iR_loc < DMR.get_nblocks(); ++iR_loc) {
//==             int iR_glob = Rgrid.mpindex.loc1D_to_glob1D(iR_loc);
//==             Peierls_phase(iR_loc) = std::exp(+im*lasA.dot(Rgrid[iR_glob]));
//==         }
//== 
//== #pragma omp parallel for schedule(static) collapse(3)
//==         for (int iblock = 0; iblock < DMR.get_nblocks(); ++iblock) {
//==             for (int irow = 0; irow < DMR.get_nrows(); ++irow) {
//==                 for (int icol = 0; icol < DMR.get_ncols(); ++icol) {
//==                     DMR(iblock, irow, icol) *= Peierls_phase(iblock);
//==                 }
//==             }
//==         }
//==        aux_DM.go_to_k();
//==    }
    


}

/// @brief Recap of all the variables of the simulation, as read from the input json file or
/// got from default values.
void Simulation::print_recap()
{
    output::title("INPUT RECAP");
    output::print("tb_model                 *", std::string(8, ' '), ctx_->cfg().tb_file());
    output::print("Space for gradient       *", std::string(8, ' '), (SpaceOfPropagation_Gradient_ == R? "R" : "k"));
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
    output::print("r0x                      *", coulomb_.get_r0()[0], " a.u.",
                                                Convert( coulomb_.get_r0()[0], AuLength, Angstrom), " angstrom");
    output::print("r0y                      *", coulomb_.get_r0()[1], " a.u.",
                                                Convert( coulomb_.get_r0()[1], AuLength, Angstrom), " angstrom");
    output::print("r0z                      *", coulomb_.get_r0()[2], " a.u.",
                                                Convert( coulomb_.get_r0()[2], AuLength, Angstrom), " angstrom");
    output::print("r0_avg                   *", coulomb_.get_r0_avg(), " a.u.",
                                                Convert( coulomb_.get_r0_avg(), AuLength, Angstrom), " angstrom");
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

// == int Simulation::get_it(const double& time__) const
// == {
// ==     return int(time__ / DEsolver_DM_.get_ResolutionTime());
// == }

/// @brief (WARNING: this function could be unsafe in the future). Get a (unique) index
/// for the creation of a node in the h5 file for the current time step.
/// @param time__ Time at which we are printing the h5 files
/// @return Index corresponding to time__ in h5 files
int Simulation::get_it_sparse(const double& time__) const
{
    static int counter = -1;
    counter++;
    return counter;
    // return int(round(time__/DEsolver_DM_.get_ResolutionTime()/ctx_->cfg().printresolution()));
}

/// @brief (For debugging and currently not used) print grids of the simulation.
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

/// @brief Check what is defined in the json file: frequency or wavelength. What is not defined has
/// the value 0.
/// @param idx__ Index of the laser for which we are checking if is defined with frequency or with wavelength
/// @return "wavelength" or "frequency", depending which one defines the laser.
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



double min(const std::vector<mdarray<double,1>>& Vec__)
{
    auto min = 1.e+07;
    for( int i = 0; i < Vec__.size(); ++i ) {
        for( int y = 0; y < Vec__[i].get_Size()[0]; ++y ) {
            min = std::min( min, *std::min(Vec__[i].begin(), Vec__[i].end()));
        }
    }
    return min;
}


double max(const std::vector<mdarray<double,1>>& Vec__)
{
    auto max = 1.e-07;
    for( int i = 0; i < Vec__.size(); ++i ) {
        for( int y = 0; y < Vec__[i].get_Size()[0]; ++y ) {
            std::vector<double> vec_i(Vec__[i].end()-Vec__[i].begin());
            max = std::max( max, *std::max_element(Vec__[i].begin(), Vec__[i].end()));
        }
    }
    return max;
}

/// @brief Print the bandstructure given a tb Hamiltonian and a kpath where we interpolate the hamiltonian
/// @param bare_kpath List of k points in lattice coordinates
/// @param Hamiltonian Hamiltonian that we diagonalize to get the band structure
void print_bandstructure(const std::vector<std::vector<double>>& bare_kpath__, Operator<std::complex<double>> Hamiltonian__)
{
    /* create the path of kpoints where to print the band structure */
    std::vector<Coordinate> path;
    path.resize(bare_kpath__.size());
    for( int ik = 0; ik < bare_kpath__.size(); ++ik ) {
        auto& bare_k = bare_kpath__[ik];
        path[ik] = Coordinate(bare_k[0], bare_k[1], bare_k[2], LatticeVectors(Space::k));
    }

    /* create kmesh */
    MeshGrid MeshGridPath(Space::k, path, 0.01);

    /* diagonalize the hamiltonian on the kmesh */
    Hamiltonian__.dft(MeshGridPath.get_mesh(), +1, false);
    std::vector<mdarray<double,1>> Eigenvalues;
    BlockMatrix<std::complex<double>> Eigenvectors;
    Hamiltonian__.get_Operator_k().diagonalize(Eigenvalues, Eigenvectors);

    /* print eigenvalues */
    std::ofstream Output;
    Output.open("BANDSTRUCTURE.txt");
    Output << "#k-number    energy(eV)\n";
    for(int ik=0; ik<Eigenvalues.size(); ik++){
        for(int iband=0; iband<Eigenvalues[ik].get_Size(0); ++iband){
            Output << std::setw(6) << ik;
            Output << std::setw(15) << std::setprecision(6) << Convert(Eigenvalues[ik](iband),AuEnergy,ElectronVolt) << std::endl;
        }
    }
    Output.close();

    std::ofstream Gnuplot;
    Gnuplot.open("plotbands.gnu");
    Gnuplot << "set xrange [ "<< -double(Eigenvalues.size())/10. << " : " << 11./10.*double(Eigenvalues.size()) << "]"<<std::endl;
    auto min_eig = Convert(min(Eigenvalues), AuEnergy, ElectronVolt);
    auto max_eig = Convert(max(Eigenvalues), AuEnergy, ElectronVolt);
    auto eig_range = max_eig - min_eig;

    Gnuplot << "set yrange [" << min_eig - eig_range/10. << ": " << max_eig + eig_range/10. << "]" << std::endl;
    int kpt;
    for( int ik = 0; ik < bare_kpath__.size(); ++ik ) {
        kpt = MeshGridPath.find(path[ik], ik);
        Gnuplot << "set arrow from  " << kpt << ", " << min_eig- eig_range/10. << " to " << kpt << ", " << max_eig+ eig_range/10. << " nohead" << std::endl;
    }
    Gnuplot << "plot \"BANDSTRUCTURE.txt\" w p" << std::endl;
    Gnuplot << "pause -1" << std::endl;
}

template<typename Func>
double findextreme(Func func, const std::vector<mdarray<double,1>>& array, int bandindex, mpi::Communicator& comm)
{
    std::vector<double> local_extreme(comm.size());
    local_extreme[kpool_comm.rank()] = findextreme(func, array, bandindex);
    MPI_Allgather(MPI_IN_PLACE, 1, MPI_DOUBLE, &local_extreme[0], 1, MPI_DOUBLE, comm.communicator());
    return *func(local_extreme.begin(), local_extreme.end());
}

template<typename Func>
double findextreme(Func func, const std::vector<mdarray<double,1>>& array, int bandindex)
{
    std::vector<double> array_k(array.size());
    for( int ik=0; ik<array.size(); ++ik ) {
        array_k[ik] = array[ik](bandindex);
    }
    return *func(array_k.begin(), array_k.end());
}

void Simulation::OpenGap()
{
    /* find maximum of valence (we know eigenvalues are sorted) and min of conduction */
    //output::print("find bangap...");
    //using fwdit = std::vector<double>::iterator;
    //auto max_valence = findextreme([](fwdit a, fwdit b){return std::max_element(a,b);},
    //                                Band_energies_, ctx_->cfg().filledbands()-1, kpool_comm);
    //auto min_conduction = findextreme([](fwdit a, fwdit b){return std::min_element(a,b);},
    //                                Band_energies_, ctx_->cfg().filledbands(), kpool_comm);
    //auto dft_bandgap = min_conduction - max_valence;
    //output::print("bandgap:", Convert(dft_bandgap, AuEnergy, ElectronVolt));
    /* open the gap with the desired value */
    //if( ctx_->cfg().opengap() < dft_bandgap && std::abs(ctx_->cfg().opengap() - dft_bandgap) > 1.e-05 ) {
    //    output::print("Warning: You want to open a gap but you are closing it!");
    //    output::print("dft band gap: ", dft_bandgap, "eV;   gap: ", ctx_->cfg().opengap(), "eV");
    //}

    Operator<std::complex<double>> Corrected_hamiltonian;
    Corrected_hamiltonian.initialize_fft(DensityMatrix_);
    Corrected_hamiltonian.lock_space(Space::k);
    Corrected_hamiltonian.lock_gauge(bloch);

    auto& Corrected_hamiltonian_k = Corrected_hamiltonian.get_Operator(Space::k);
    auto& Corrected_hamiltonian_R = Corrected_hamiltonian.get_Operator(Space::R);
    Corrected_hamiltonian_k.fill(0.);
    for( int ik = 0; ik < Corrected_hamiltonian_k.get_nblocks(); ++ik ) {
        for ( int ival = 0; ival < ctx_->cfg().filledbands(); ++ival ) {
            Corrected_hamiltonian_k(ik, ival, ival) = Band_energies_[ik](ival) - ctx_->cfg().opengap()/2.;
        }
        for ( int icond = ctx_->cfg().filledbands(); icond < Corrected_hamiltonian_k.get_nrows(); ++icond ) {
            Corrected_hamiltonian_k(ik, icond, icond) = Band_energies_[ik](icond) + ctx_->cfg().opengap()/2.;
        }
    }

    Corrected_hamiltonian.go_to_wannier();
    Corrected_hamiltonian.go_to_R();

    /* copy in the material hamiltonian the corrected one */
    material_.H.initialized_dft = false;
    material_.H.initialize_fft(DensityMatrix_);
    std::copy(Corrected_hamiltonian_k.begin(), Corrected_hamiltonian_k.end(), material_.H.get_Operator_k().begin());
    std::copy(Corrected_hamiltonian_R.begin(), Corrected_hamiltonian_R.end(), material_.H.get_Operator_R().begin());

}


void Simulation::Apply_Peierls_phase(Operator<std::complex<double>>& O__, const double& time__, const int sign = +1)
{
    O__.go_to_R();
    static mdarray<std::complex<double>,1> Peierls_phase({O__.get_Operator(Space::R).get_nblocks()});

    /* Calculate Peierls phase on the grid */
    auto At     = setoflaser_.VectorPotential(time__);
    auto& Rgrid = MeshGrid::MasterRgrid_GammaCentered; //WARNING! Here we are supposing O__ R grid is the MasterRgrid! A check would be ideal

#pragma omp parallel for schedule(static)
    for (int iR_loc = 0; iR_loc < O__.get_Operator(Space::R).get_nblocks(); ++iR_loc) {
        int iR_glob = Rgrid.mpindex.loc1D_to_glob1D(iR_loc);
        Peierls_phase(iR_loc) = std::exp(im*double(sign)*At.dot(Rgrid[iR_glob]));
    }

#pragma omp parallel for schedule(static) collapse(3)
    for (int iblock = 0; iblock < O__.get_Operator(Space::R).get_nblocks(); ++iblock) {
        for (int irow = 0; irow < O__.get_Operator(Space::R).get_nrows(); ++irow) {
            for (int icol = 0; icol < O__.get_Operator(Space::R).get_ncols(); ++icol) {
                O__.get_Operator(Space::R)(iblock, irow, icol) *= Peierls_phase(iblock);
            }
        }
    }
}
