#include <cmath>
#include <filesystem>
#include <limits>
#include <set>
#include <sstream>
#include "Simulation/Simulation.hpp"
#include "core/mpi/Communicator.hpp"
#ifdef EDUS_HDF5
#include "core/hdf5/hdf5_tree.hpp"
#endif

/* Checkpoint of the time propagation, to restart a simulation from a given step.
   The state at the beginning of the step n (before printing the outputs at t_n) is given by:
   - the density matrix rho(k) in the wannier gauge (only the k part is needed: the R part is
     recomputed from it in the source term)
   - for AB of order p: the source terms f(t_{n-i}, rho(t_{n-i})), i=1..p-1
   - for each laser: the vector potential A and the time of the last call of VectorPotential,
     since A is integrated step by step (Laser::VectorPotential)
   - the variables of the outputs that depend on the previous steps (work of the field, h5 index, ...)
   H0, DM0, the interaction and the lasers' fields are recomputed from the input.
   The k arrays are written with global k indices: the restart can use a different number of MPI ranks. */

namespace
{
    constexpr int checkpoint_version = 1;
    constexpr double nan_ = std::numeric_limits<double>::quiet_NaN();

    /// Keys of the input that can change between the run that wrote the checkpoint and the restart
    const std::set<std::string> free_keys = {
        "lasers", "finaltime", "finaltime_units", "printresolution", "printresolution_pulse",
        "printresolution_h5", "toprint", "kpath", "restart_file", "restart_check_input",
        "checkpoint_resolution", "checkpoint_times", "checkpoint_times_units", "locked" };

    /// Parameters that define a laser: t0, omega, amplitude, cycles, phase, polarization
    std::vector<double> laser_signature(Laser& laser__)
    {
        auto& pol = laser__.get_Polarization().get("Cartesian");
        return { laser__.get_InitialTime(), laser__.get_Omega(), laser__.get_Amplitude(),
                 laser__.get_NumberOfCycles(), laser__.get_Phase(), pol[0], pol[1], pol[2] };
    }

    bool same_signature(const std::vector<double>& a__, const std::vector<double>& b__)
    {
        for ( size_t i = 0; i < a__.size(); ++i ) {
            if ( std::abs(a__[i] - b__[i]) > 1.e-12 * std::max(1., std::abs(a__[i])) ) {
                return false;
            }
        }
        return true;
    }

    bool is_root() { return mpi::Communicator::world().rank() == 0; }

    template <typename T>
    void bcast(T* data__, int count__)
    {
#ifdef EDUS_MPI
        MPI_Bcast(data__, count__, mpi::type_wrapper<T>(), 0, mpi::Communicator::world().communicator());
#endif
    }

    void bcast(std::string& str__)
    {
        int size = str__.size();
        bcast(&size, 1);
        str__.resize(size);
        bcast(&str__[0], size);
    }

    /// Copies the local k blocks of O__ in global__ (k global index, row, col), allocated only on rank 0
    void gather_k(const BlockMatrix<std::complex<double>>& O__, const MPIindex<3>& mpindex__, int num_kpoints__,
                  std::vector<std::complex<double>>& global__)
    {
        int nrows = O__.get_nrows();
        int ncols = O__.get_ncols();
        int nlocal = mpindex__.get_nlocal();
        std::vector<std::complex<double>> local(size_t(nlocal) * nrows * ncols);
        #pragma omp parallel for
        for ( int ik = 0; ik < nlocal; ++ik ) {
            for ( int irow = 0; irow < nrows; ++irow ) {
                for ( int icol = 0; icol < ncols; ++icol ) {
                    local[(size_t(ik) * nrows + irow) * ncols + icol] = O__(ik, irow, icol);
                }
            }
        }
#ifdef EDUS_MPI
        auto& comm = mpi::Communicator::world();
        int count = nlocal * nrows * ncols;
        int offset = int(mpindex__.get_LocalRange().first) * nrows * ncols;
        std::vector<int> counts(comm.size()), offsets(comm.size());
        MPI_Gather(&count, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, comm.communicator());
        MPI_Gather(&offset, 1, MPI_INT, offsets.data(), 1, MPI_INT, 0, comm.communicator());
        global__.resize( is_root() ? size_t(num_kpoints__) * nrows * ncols : 0 );
        MPI_Gatherv(local.data(), count, MPI_C_DOUBLE_COMPLEX, global__.data(), counts.data(), offsets.data(),
                    MPI_C_DOUBLE_COMPLEX, 0, comm.communicator());
#else
        global__ = local;
#endif
    }

    /// Inverse of gather_k: copies in the local k blocks of O__ the values of global__ (given on rank 0)
    void scatter_k(const std::vector<std::complex<double>>& global__, const MPIindex<3>& mpindex__,
                   BlockMatrix<std::complex<double>>& O__)
    {
        int nrows = O__.get_nrows();
        int ncols = O__.get_ncols();
        int nlocal = mpindex__.get_nlocal();
        std::vector<std::complex<double>> local(size_t(nlocal) * nrows * ncols);
#ifdef EDUS_MPI
        auto& comm = mpi::Communicator::world();
        int count = nlocal * nrows * ncols;
        int offset = int(mpindex__.get_LocalRange().first) * nrows * ncols;
        std::vector<int> counts(comm.size()), offsets(comm.size());
        MPI_Gather(&count, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, comm.communicator());
        MPI_Gather(&offset, 1, MPI_INT, offsets.data(), 1, MPI_INT, 0, comm.communicator());
        MPI_Scatterv(global__.data(), counts.data(), offsets.data(), MPI_C_DOUBLE_COMPLEX,
                     local.data(), count, MPI_C_DOUBLE_COMPLEX, 0, comm.communicator());
#else
        local = global__;
#endif
        #pragma omp parallel for
        for ( int ik = 0; ik < nlocal; ++ik ) {
            for ( int irow = 0; irow < nrows; ++irow ) {
                for ( int icol = 0; icol < ncols; ++icol ) {
                    O__(ik, irow, icol) = local[(size_t(ik) * nrows + irow) * ncols + icol];
                }
            }
        }
    }

    [[noreturn]] void checkpoint_error(const std::string& name__, const std::string& msg__)
    {
        throw std::runtime_error("restart from " + name__ + ": " + msg__);
    }
}

/// @brief Writes the checkpoint of the current step in name__ (through a temporary file, so that a
/// previous checkpoint with the same name is replaced only when the new one is complete).
/// To be called at the beginning of a step, before printing the outputs.
void Simulation::write_checkpoint(const std::string& name__)
{
#ifndef EDUS_HDF5
    throw std::runtime_error("checkpoints need EDUS compiled with HDF5 (EDUS_HDF5=ON)");
#else
    PROFILE("Simulation::write_checkpoint");
    if ( propagator_.processor() == device ) {
        throw std::runtime_error("checkpoints are not implemented on gpu");
    }
    auto& desolver = propagator_.desolver();
    if ( !desolver.history_complete() ) {
        output::print("checkpoint not written at step ", std::to_string(desolver.istep()), ": AB history not complete yet");
        return;
    }
    auto& mpindex = decomposition_.mpindex();
    auto& DM = electrons_state_.DensityMatrix();
    int num_bands = DM.get_Operator_k().get_nrows();
    int num_kpoints = gridstructure_.kgrid()->get_TotalSize();
    std::string tmp = name__ + ".tmp";

    std::unique_ptr<HDF5_tree> fout;
    if ( is_root() ) {
        fout = std::make_unique<HDF5_tree>(tmp, hdf5_access_t::truncate);
        fout->write("version", checkpoint_version);
        fout->write("istep", double(desolver.istep()));
        fout->write("time", desolver.current_time());
        fout->write("dt", desolver.time_step());
        fout->write("initial_time", desolver.parameters().initial_time);
        fout->write("solver", int(desolver.parameters().solver));
        fout->write("order", desolver.parameters().order);
        auto kgrid = ctx_->cfg().grid();
        fout->write("grid", std::vector<int>(kgrid.begin(), kgrid.end()));
        fout->write("num_bands", num_bands);
        fout->write("num_kpoints", num_kpoints);
        fout->write("num_history", desolver.num_saved_history());
        fout->write("input", ctx_->cfg().dict().dump());

        auto s = output_.restart_state();
        fout->write("output_state", std::vector<double>{ s.work, s.last_power, s.last_time, s.initial_energy,
                                                         double(s.first_energy_step), double(s.index_h5),
                                                         double(s.num_prints) });
        fout->write("num_lasers", setoflaser_.size());
        fout->create_node("lasers");
        for ( int ilaser = 0; ilaser < setoflaser_.size(); ++ilaser ) {
            auto& laser = setoflaser_[ilaser];
            auto& A = laser.get_A().get("Cartesian");
            (*fout)["lasers"].create_node(ilaser);
            (*fout)["lasers"][ilaser].write("signature", laser_signature(laser));
            (*fout)["lasers"][ilaser].write("A_state", std::vector<double>{ A[0], A[1], A[2], laser.get_previous_Time() });
        }
    }

    /* k arrays, one at a time to limit the memory on rank 0 */
    std::vector<std::complex<double>> global;
    auto write_k = [&](const BlockMatrix<std::complex<double>>& O, const std::string& label) {
        gather_k(O, mpindex, num_kpoints, global);
        if ( is_root() ) {
            fout->write(label, reinterpret_cast<double*>(global.data()), int(2 * global.size()));
        }
    };
    write_k(DM.get_Operator_k(), "DMk");
    for ( int i = 1; i <= desolver.num_saved_history(); ++i ) {
        write_k(desolver.history(i).get_Operator_k(), "history_" + std::to_string(i));
    }

    if ( is_root() ) {
        fout.reset();
        std::filesystem::rename(tmp, name__);
    }
    mpi::Communicator::world().barrier();
    output::print("-> checkpoint written in ", name__, " (step ", std::to_string(desolver.istep()), ", time ",
                  Convert(desolver.current_time(), AuTime, FemtoSeconds), " fs)");
#endif
}

/// @brief Reads the step of the checkpoint and the number of steps printed in the txt files before it,
/// needed to initialize the outputs before the rest of the restart
void Simulation::read_checkpoint_header(const std::string& name__, long& istep__, long& num_prints__)
{
#ifndef EDUS_HDF5
    throw std::runtime_error("restart needs EDUS compiled with HDF5 (EDUS_HDF5=ON)");
#else
    int ok = 1;
    double header[2] = { nan_, nan_ };
    if ( is_root() ) {
        ok = std::filesystem::exists(name__);
        if ( ok ) {
            HDF5_tree fin(name__, hdf5_access_t::read_only);
            int version = -1;
            fin.read("version", &version, 1);
            std::vector<double> s(7, nan_);
            fin.read("output_state", s);
            fin.read("istep", &header[0], 1);
            header[1] = s[6];
            ok = ( version == checkpoint_version );
        }
    }
    bcast(&ok, 1);
    bcast(header, 2);
    if ( !ok || std::isnan(header[0]) || std::isnan(header[1]) ) {
        checkpoint_error(name__, "file not found or not a checkpoint of this version of EDUS");
    }
    istep__ = std::lround(header[0]);
    num_prints__ = std::lround(header[1]);
#endif
}

/// @brief Sets the density matrix, the solver, the lasers and the outputs to the state of the checkpoint.
/// Lasers: each laser of the input that starts before the restart time must be in the checkpoint (same
/// parameters): its vector potential is restored. Lasers starting after the restart time are new.
void Simulation::restore_checkpoint(const std::string& name__)
{
#ifndef EDUS_HDF5
    throw std::runtime_error("restart needs EDUS compiled with HDF5 (EDUS_HDF5=ON)");
#else
    PROFILE("Simulation::restore_checkpoint");
    if ( propagator_.processor() == device ) {
        throw std::runtime_error("restart is not implemented on gpu");
    }
    auto& desolver = propagator_.desolver();
    auto& mpindex = decomposition_.mpindex();
    auto& DM = electrons_state_.DensityMatrix();
    int num_bands = DM.get_Operator_k().get_nrows();
    int num_kpoints = gridstructure_.kgrid()->get_TotalSize();

    /* scalars: read on rank 0, broadcast, checked on all ranks (so that all of them throw) */
    std::vector<double> scalars(4, nan_);   /* istep, time, dt, initial_time */
    std::vector<int> iscalars(6, -1);      /* solver, order, num_bands, num_kpoints, num_history, num_lasers */
    std::vector<int> grid(3, -1);
    std::vector<double> out_state(7, nan_);
    std::string input;
    std::vector<double> lasers_data;       /* 12 values per laser: signature (8), A (3), previous time */
    std::unique_ptr<HDF5_tree> fin;
    if ( is_root() ) {
        fin = std::make_unique<HDF5_tree>(name__, hdf5_access_t::read_only);
        fin->read("istep", &scalars[0], 1);
        fin->read("time", &scalars[1], 1);
        fin->read("dt", &scalars[2], 1);
        fin->read("initial_time", &scalars[3], 1);
        fin->read("solver", &iscalars[0], 1);
        fin->read("order", &iscalars[1], 1);
        fin->read("num_bands", &iscalars[2], 1);
        fin->read("num_kpoints", &iscalars[3], 1);
        fin->read("num_history", &iscalars[4], 1);
        fin->read("num_lasers", &iscalars[5], 1);
        fin->read("grid", grid);
        fin->read("output_state", out_state);
        fin->read("input", input);
        lasers_data.assign(12 * std::max(iscalars[5], 0), nan_);
        for ( int ilaser = 0; ilaser < iscalars[5]; ++ilaser ) {
            std::vector<double> sig(8, nan_), st(4, nan_);
            (*fin)["lasers"][ilaser].read("signature", sig);
            (*fin)["lasers"][ilaser].read("A_state", st);
            std::copy(sig.begin(), sig.end(), lasers_data.begin() + 12 * ilaser);
            std::copy(st.begin(), st.end(), lasers_data.begin() + 12 * ilaser + 8);
        }
    }
    bcast(scalars.data(), 4);
    bcast(iscalars.data(), 6);
    bcast(grid.data(), 3);
    bcast(out_state.data(), 7);
    bcast(input);
    int nlaser_data = lasers_data.size();
    bcast(&nlaser_data, 1);
    lasers_data.resize(nlaser_data);
    bcast(lasers_data.data(), nlaser_data);

    for ( auto x : scalars ) {
        if ( std::isnan(x) ) checkpoint_error(name__, "missing data in the file");
    }
    for ( auto x : out_state ) {
        if ( std::isnan(x) ) checkpoint_error(name__, "missing output state in the file");
    }
    for ( auto x : lasers_data ) {
        if ( std::isnan(x) ) checkpoint_error(name__, "missing laser data in the file");
    }

    /* the propagation must be the same */
    long istep = std::lround(scalars[0]);
    auto& dp = desolver.parameters();
    std::stringstream err;
    if ( std::abs(scalars[2] - dp.dt) > 1.e-14 * dp.dt ) err << "dt " << scalars[2] << " vs " << dp.dt << " a.u.; ";
    if ( std::abs(scalars[3] - dp.initial_time) > 1.e-12 * std::max(1., std::abs(dp.initial_time)) ) {
        err << "initial time " << scalars[3] << " vs " << dp.initial_time << " a.u.; ";
    }
    if ( iscalars[0] != int(dp.solver) || iscalars[1] != dp.order ) err << "solver or order; ";
    if ( iscalars[2] != num_bands ) err << "number of bands " << iscalars[2] << " vs " << num_bands << "; ";
    if ( iscalars[3] != num_kpoints ) err << "number of k points " << iscalars[3] << " vs " << num_kpoints << "; ";
    auto kgrid = ctx_->cfg().grid();
    if ( grid[0] != kgrid[0] || grid[1] != kgrid[1] || grid[2] != kgrid[2] ) err << "k grid; ";
    if ( iscalars[4] != desolver.num_saved_history() ) err << "number of history arrays; ";
    if ( !err.str().empty() ) {
        checkpoint_error(name__, "different propagation: " + err.str());
    }

    /* the rest of the input must be the same, apart from lasers, final time and outputs */
    auto parent = nlohmann::json::parse(input);
    auto& current = ctx_->cfg().dict();
    std::stringstream diff;
    for ( auto* d : { &parent, &current } ) {
        for ( auto it = d->begin(); it != d->end(); ++it ) {
            if ( free_keys.count(it.key()) ) continue;
            auto& other = ( d == &parent ? current : parent );
            bool equal = other.contains(it.key()) && other.at(it.key()) == it.value();
            if ( !equal && diff.str().find(" " + it.key() + " ") == std::string::npos ) {
                diff << " " << it.key() << " ";
            }
        }
    }
    if ( !diff.str().empty() ) {
        if ( ctx_->cfg().restart_check_input() ) {
            checkpoint_error(name__, "input different from the one of the checkpoint for the keys:" + diff.str() +
                             "(set restart_check_input to false to restart anyway)");
        }
        output::print("WARNING: restart with input different from the checkpoint for the keys:", diff.str());
    }

    /* lasers */
    double restart_time = scalars[1];
    int num_parent_lasers = iscalars[5];
    std::vector<bool> used(num_parent_lasers, false);
    /* all the lasers are called together: a new laser starts from the time of the last call */
    double previous_time = ( num_parent_lasers > 0 ? lasers_data[11] : restart_time );
    for ( int ilaser = 0; ilaser < setoflaser_.size(); ++ilaser ) {
        auto& laser = setoflaser_[ilaser];
        auto signature = laser_signature(laser);
        int found = -1;
        for ( int jlaser = 0; jlaser < num_parent_lasers && found < 0; ++jlaser ) {
            std::vector<double> sig(lasers_data.begin() + 12 * jlaser, lasers_data.begin() + 12 * jlaser + 8);
            if ( !used[jlaser] && same_signature(signature, sig) ) found = jlaser;
        }
        if ( found >= 0 ) {
            used[found] = true;
            auto* st = &lasers_data[12 * found + 8];
            laser.set_A_state(Coordinate(st[0], st[1], st[2]), st[3]);
        }
        else if ( laser.get_InitialTime() < restart_time - 1.e-07 ) {
            checkpoint_error(name__, "laser " + std::to_string(ilaser) + " starts before the restart time but " +
                             "is not in the checkpoint");
        }
        else {
            laser.set_A_state(Coordinate(0., 0., 0.), previous_time);
        }
    }
    for ( int jlaser = 0; jlaser < num_parent_lasers; ++jlaser ) {
        if ( !used[jlaser] && lasers_data[12 * jlaser] < restart_time - 1.e-07 ) {
            checkpoint_error(name__, "laser " + std::to_string(jlaser) + " of the checkpoint started before the " +
                             "restart time and is not in the input");
        }
    }

    /* density matrix and history of the solver */
    desolver.restart_from_step(istep);
    std::vector<std::complex<double>> global;
    auto read_k = [&](BlockMatrix<std::complex<double>>& O, const std::string& label) {
        int ok = 1;
        if ( is_root() ) {
            global.assign(size_t(num_kpoints) * num_bands * num_bands, std::complex<double>(nan_, nan_));
            fin->read(label, reinterpret_cast<double*>(global.data()), int(2 * global.size()));
            ok = !std::isnan(global.back().real()) && !std::isnan(global.front().real());
        }
        bcast(&ok, 1);
        if ( !ok ) checkpoint_error(name__, "missing " + label);
        scatter_k(global, mpindex, O);
    };
    read_k(DM.get_Operator_k(), "DMk");
    DM.lock_gauge(wannier);
    DM.lock_space(k);
    for ( int i = 1; i <= desolver.num_saved_history(); ++i ) {
        auto& h = desolver.history(i);
        read_k(h.get_Operator_k(), "history_" + std::to_string(i));
        h.lock_space(k);
    }

    output_.set_restart_state({ out_state[0], out_state[1], out_state[2], out_state[3], out_state[4] > 0.5,
                                int(std::lround(out_state[5])), std::lround(out_state[6]) });
    fin.reset();
    mpi::Communicator::world().barrier();
    output::print("-> restart from ", name__, " at step ", std::to_string(istep), ", time ", Convert(restart_time, AuTime, FemtoSeconds),
                  " fs");
#endif
}
