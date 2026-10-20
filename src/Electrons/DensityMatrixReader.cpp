#include <cstdio>
#include <fstream>
#include <iomanip>
#include <set>
#include <filesystem>
#include "Electrons/DensityMatrixReader.hpp"
#include "ConvertUnits.hpp"
#include "GlobalFunctions.hpp"

namespace electron
{
    /// Rank 0 found an error: all the ranks throw, with the message of rank 0
    static void check_rank0(const std::string& error__, const mpi::Communicator& comm__)
    {
        int failed = !error__.empty();
#ifdef EDUS_MPI
        MPI_Bcast(&failed, 1, MPI_INT, 0, comm__.communicator());
#endif
        if ( failed ) {
            throw std::runtime_error(comm__.rank() == 0 ? error__ : "DensityMatrixReader: error on rank 0");
        }
    }

    std::string DensityMatrixReader::file(int i1__, int i2__, int isnap__) const
    {
        return directory_ + "/" + basename_[i1__][i2__] + "_" + std::to_string(isnap__) + ".txt";
    }

    void DensityMatrixReader::initialize(const std::string& directory__, const std::string& kpoints_file__,
                                         const std::string& times_file__, int num_wann__, const MeshGrid& kgrid__,
                                         const mpi::Communicator& comm__)
    {
        PROFILE("DensityMatrixReader::initialize");
        directory_ = directory__;
        num_wann_ = num_wann__;
        int num_k = kgrid__.get_TotalSize();

        /* file names: indices are not separated, so we check that they are unique */
        basename_.assign(num_wann_, std::vector<std::string>(num_wann_));
        std::set<std::string> names;
        for (int i1 = 0; i1 < num_wann_; ++i1) {
            for (int i2 = i1; i2 < num_wann_; ++i2) {
                basename_[i1][i2] = "P_in_W_" + std::to_string(i1) + std::to_string(i2);
                if ( !names.insert(basename_[i1][i2]).second ) {
                    throw std::runtime_error("DensityMatrixReader: file name " + basename_[i1][i2] + " is ambiguous");
                }
            }
        }

        /* map the k points of the files to the EDUS grid */
        edus_index_.assign(num_k, -1);
        std::string error;
        if ( comm__.rank() == 0 ) {
            std::ifstream is(directory_ + "/" + kpoints_file__);
            std::vector<int> found(num_k, 0);
            int iline = 0;
            double k0, k1, k2;
            if ( !is ) {
                error = "DensityMatrixReader: cannot open " + directory_ + "/" + kpoints_file__;
            }
            while ( error.empty() && is >> k0 >> k1 >> k2 ) {
                if ( iline >= num_k ) {
                    error = "DensityMatrixReader: more k points in " + kpoints_file__ + " than in the EDUS grid";
                    break;
                }
                Coordinate kfile(k0, k1, k2, LatticeVectors(k));
                int index = kgrid__.find(kfile);
                /* find rounds to the closest grid point: check that the k point is really there */
                bool on_grid = ( index >= 0 && index < num_k );
                if ( on_grid ) {
                    auto diff = ( kgrid__[index] - kfile ).get(LatticeVectors(k));
                    for (auto ix : {0, 1, 2}) {
                        on_grid = on_grid && std::abs(diff[ix] - std::round(diff[ix])) < 1.e-05;
                    }
                }
                if ( !on_grid || found[index] ) {
                    std::stringstream ss;
                    ss << "DensityMatrixReader: k point " << iline << " (" << k0 << ", " << k1 << ", " << k2 << ") "
                       << ( on_grid ? "appears twice" : "is not on the EDUS k grid" );
                    error = ss.str();
                    break;
                }
                found[index] = 1;
                edus_index_[iline] = index;
                ++iline;
            }
            if ( error.empty() && iline != num_k ) {
                error = "DensityMatrixReader: " + std::to_string(iline) + " k points in " + kpoints_file__ +
                        ", while the EDUS grid has " + std::to_string(num_k);
            }
        }
        check_rank0(error, comm__);
#ifdef EDUS_MPI
        MPI_Bcast(edus_index_.data(), num_k, MPI_INT, 0, comm__.communicator());
#endif

        /* snapshots: P_in_W_00_0.txt, P_in_W_00_1.txt, ... */
        if ( comm__.rank() == 0 ) {
            num_snapshots_ = 0;
            while ( std::filesystem::exists(file(0, 0, num_snapshots_)) ) {
                ++num_snapshots_;
            }
            if ( num_snapshots_ == 0 ) {
                error = "DensityMatrixReader: " + file(0, 0, 0) + " not found";
            }
        }
        check_rank0(error, comm__);
#ifdef EDUS_MPI
        MPI_Bcast(&num_snapshots_, 1, MPI_INT, 0, comm__.communicator());
#endif

        /* times of the snapshots (optional) */
        times_.assign(num_snapshots_, 0.);
        int has_times = 0;
        if ( comm__.rank() == 0 && !times_file__.empty() ) {
            std::ifstream is(directory_ + "/" + times_file__);
            int it = 0;
            while ( it < num_snapshots_ && is >> times_[it] ) {
                ++it;
            }
            has_times = ( it == num_snapshots_ );
        }
#ifdef EDUS_MPI
        MPI_Bcast(&has_times, 1, MPI_INT, 0, comm__.communicator());
        MPI_Bcast(times_.data(), num_snapshots_, MPI_DOUBLE, 0, comm__.communicator());
#endif
        if ( !has_times ) {
            for (int it = 0; it < num_snapshots_; ++it) {
                times_[it] = it;
            }
        }
        output::print("-> density matrix from files in ", directory_, ": ", num_snapshots_, " snapshots, ", num_k, " k points");
        if ( !has_times ) {
            output::print("   no times found in ", times_file__, ", snapshots labelled with their index");
        }
    }

    void DensityMatrixReader::read(int isnap__, Operator<std::complex<double>>& DM__, const MPIindex<3>& mpindex__,
                                   const mpi::Communicator& comm__) const
    {
        PROFILE("DensityMatrixReader::read");
        int num_k = int(edus_index_.size());
        std::vector<std::complex<double>> values(num_k);
        auto& DMk = DM__.get_Operator(Space::k);
        int nlocal = DMk.get_nblocks();

        for (int i1 = 0; i1 < num_wann_; ++i1) {
            for (int i2 = i1; i2 < num_wann_; ++i2) {
                std::string error;
                if ( comm__.rank() == 0 ) {
                    auto filename = file(i1, i2, isnap__);
                    std::FILE* fp = std::fopen(filename.c_str(), "r");
                    if ( fp == nullptr ) {
                        error = "DensityMatrixReader: cannot open " + filename;
                    }
                    else {
                        double re, im;
                        for (int iline = 0; iline < num_k; ++iline) {
                            if ( std::fscanf(fp, "%lf %lf", &re, &im) != 2 ) {
                                error = "DensityMatrixReader: " + filename + " has less than " + std::to_string(num_k) + " lines";
                                break;
                            }
                            values[edus_index_[iline]] = std::complex<double>(re, im);
                        }
                        std::fclose(fp);
                    }
                }
                check_rank0(error, comm__);
#ifdef EDUS_MPI
                MPI_Bcast(values.data(), num_k, MPI_CXX_DOUBLE_COMPLEX, 0, comm__.communicator());
#endif
                /* the file contains rho_{i2,i1}: the other triangle from hermiticity */
                for (int ik = 0; ik < nlocal; ++ik) {
                    auto& value = values[mpindex__.loc1D_to_glob1D(ik)];
                    DMk(ik, i2, i1) = value;
                    DMk(ik, i1, i2) = std::conj(value);
                }
            }
        }
        DM__.lock_gauge(wannier);
        DM__.lock_space(Space::k);
    }

    /// Gathers on rank 0 a block matrix distributed in contiguous slabs of blocks (k or R) among the ranks
    static std::vector<std::complex<double>> gather_blocks(const BlockMatrix<std::complex<double>>& O__,
                                                           const mpi::Communicator& comm__)
    {
        int count = O__.get_nblocks() * O__.get_nrows() * O__.get_ncols();
        std::vector<int> counts(comm__.size(), count), displs(comm__.size(), 0);
#ifdef EDUS_MPI
        MPI_Gather(&count, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, comm__.communicator());
#endif
        int total = 0;
        for (int irank = 0; irank < comm__.size(); ++irank) {
            displs[irank] = total;
            total += counts[irank];
        }
        std::vector<std::complex<double>> global(comm__.rank() == 0 ? total : 0);
#ifdef EDUS_MPI
        MPI_Gatherv(&(O__(0,0,0)), count, MPI_CXX_DOUBLE_COMPLEX,
                    global.data(), counts.data(), displs.data(), MPI_CXX_DOUBLE_COMPLEX, 0, comm__.communicator());
#else
        std::copy(&(O__(0,0,0)), &(O__(0,0,0)) + count, global.begin());
#endif
        return global;
    }

    void print_population_k(const std::string& filename__, const BlockMatrix<std::complex<double>>& DMk_bloch__,
                            const BandStructure& bandstructure__, const MeshGrid& kgrid__,
                            const MPIindex<3>& mpindex__, int filledbands__, const mpi::Communicator& comm__)
    {
        PROFILE("print_population_k");
        int nbnd = DMk_bloch__.get_nrows();
        int nlocal = DMk_bloch__.get_nblocks();

        /* pack populations and energies of the local k points, to gather them with the same routine */
        BlockMatrix<std::complex<double>> local(k, nlocal, 2, nbnd);
        for (int ik = 0; ik < nlocal; ++ik) {
            for (int ib = 0; ib < nbnd; ++ib) {
                local(ik, 0, ib) = DMk_bloch__(ik, ib, ib);
                local(ik, 1, ib) = bandstructure__.energies()[ik](ib);
            }
        }
        auto global = gather_blocks(local, comm__);
        if ( comm__.rank() != 0 ) {
            return;
        }
        int num_k = kgrid__.get_TotalSize();
        double min_pop = 1., max_pop = 0., max_imag = 0.;
        std::ofstream os(filename__);
        os << "# populations in the bloch gauge P_n(k) = rho_nn(k), on the EDUS k grid\n";
        os << "# k1 k2 k3 (crystal) | E_n (eV), n = 1.." << nbnd << " | P_n, n = 1.." << nbnd << "\n";
        for (int ik = 0; ik < num_k; ++ik) {
            auto& kcrys = kgrid__[ik].get(LatticeVectors(k));
            os << std::setw(12) << std::setprecision(6) << kcrys[0] << std::setw(12) << kcrys[1] << std::setw(12) << kcrys[2];
            for (int ib = 0; ib < nbnd; ++ib) {
                os << std::setw(14) << std::setprecision(6) << Convert(global[(ik * 2 + 1) * nbnd + ib].real(), AuEnergy, ElectronVolt);
            }
            for (int ib = 0; ib < nbnd; ++ib) {
                auto pop = global[ik * 2 * nbnd + ib];
                os << std::setw(22) << std::setprecision(12) << pop.real();
                min_pop = std::min(min_pop, pop.real());
                max_pop = std::max(max_pop, pop.real());
                max_imag = std::max(max_imag, std::abs(pop.imag()));
            }
            os << "\n";
        }
        output::print("   P_n(k) in [", min_pop, ", ", max_pop, "], max |Im P_n(k)| = ", max_imag, " -> ", filename__);
    }

    void print_fatbands_population(const std::string& filename__, const BlockMatrix<std::complex<double>>& DMR__,
                                   Operator<std::complex<double>> H__, const std::vector<std::vector<double>>& kpath__,
                                   const GridStructure& gridstructure__, int filledbands__, double opengap__,
                                   const mpi::Communicator& comm__)
    {
        PROFILE("print_fatbands_population");
        int nbnd = DMR__.get_nrows();
        auto global = gather_blocks(DMR__, comm__);
        if ( comm__.rank() != 0 ) {
            return;
        }

        /* rho(R) on the gamma centered R grid. For even sizes, R_i = -N_i/2 and +N_i/2 are the same point of the
           grid: we split it in the two images with half weight, so that rho(-R) = rho(R)^dagger and the
           interpolated rho(k) is hermitian */
        auto& Rgrid = *gridstructure__.Rgrid_GammaCentered();
        auto& size = gridstructure__.size();
        int nR = Rgrid.get_TotalSize();
        std::vector<std::array<double,3>> Rvectors;
        std::vector<int> Rindex;
        std::vector<double> weight;
        for (int iR = 0; iR < nR; ++iR) {
            auto Rcrys = Rgrid[iR].get(LatticeVectors(R));
            std::vector<std::array<double,3>> images = { {Rcrys[0], Rcrys[1], Rcrys[2]} };
            for (auto ix : {0, 1, 2}) {
                if ( size[ix] % 2 == 0 && std::abs(std::abs(Rcrys[ix]) - size[ix] / 2.) < 1.e-06 ) {
                    int nimages = int(images.size());
                    for (int im = 0; im < nimages; ++im) {
                        auto image = images[im];
                        image[ix] = -image[ix];
                        images.push_back(image);
                    }
                }
            }
            for (auto& image : images) {
                Rvectors.push_back(image);
                Rindex.push_back(iR);
                weight.push_back(1. / images.size());
            }
        }
        int nR_ext = int(Rvectors.size());
        mdarray<std::complex<double>,3> values({nR_ext, nbnd, nbnd});
        mdarray<double,2> bare_mesh({nR_ext, 3});
        for (int iR = 0; iR < nR_ext; ++iR) {
            for (int i1 = 0; i1 < nbnd; ++i1) {
                for (int i2 = 0; i2 < nbnd; ++i2) {
                    values(iR, i1, i2) = weight[iR] * global[(Rindex[iR] * nbnd + i1) * nbnd + i2];
                }
            }
            for (auto ix : {0, 1, 2}) {
                bare_mesh(iR, ix) = Rvectors[iR][ix];
            }
        }
        Operator<std::complex<double>> rho;
        rho.get_Operator_R().initialize(R, values);
        rho.get_Operator_R().set_MeshGrid(MeshGrid(R, bare_mesh, LatticeVectors(R)));
        rho.lock_gauge(wannier);
        rho.lock_space(R);

        /* path, as in print_bandstructure_path */
        std::vector<Coordinate> vertices(kpath__.size());
        for (int ivertex = 0; ivertex < int(kpath__.size()); ++ivertex) {
            vertices[ivertex] = Coordinate(kpath__[ivertex][0], kpath__[ivertex][1], kpath__[ivertex][2], LatticeVectors(Space::k));
        }
        MeshGrid MeshGridPath(Space::k, vertices, 0.01);

        /* interpolation on the path: bands of H0 and density matrix in the wannier gauge */
        H__.dft(MeshGridPath.get_mesh(), +1, false);
        BandStructure bandstructure;
        bandstructure.initialize(H__.get_Operator_k());
        std::vector<double> shift(nbnd);
        for (int ib = 0; ib < nbnd; ++ib) {
            shift[ib] = ( ib < filledbands__ ? -1. : +1. ) * opengap__ / 2.;
        }
        bandstructure.shift_energies(shift);
        auto& energies = bandstructure.energies();
        auto& U = bandstructure.U();
        rho.dft(MeshGridPath.get_mesh(), +1, false);
        auto& rhok = rho.get_Operator_k();

        /* P_n(k) = (U^dagger rho U)_nn; carriers: holes (1-P) in filled bands, electrons (P) in empty ones */
        int num_k = int(energies.size());
        std::vector<double> population(num_k * nbnd), carriers(num_k * nbnd);
        double max_carriers = 0., min_eig = energies[0](0), max_eig = energies[0](0);
        for (int ik = 0; ik < num_k; ++ik) {
            for (int ib = 0; ib < nbnd; ++ib) {
                std::complex<double> value = 0.;
                for (int i1 = 0; i1 < nbnd; ++i1) {
                    for (int i2 = 0; i2 < nbnd; ++i2) {
                        value += std::conj(U[ik](i1, ib)) * rhok(ik, i1, i2) * U[ik](i2, ib);
                    }
                }
                population[ik * nbnd + ib] = value.real();
                carriers[ik * nbnd + ib] = ( ib < filledbands__ ? 1. - value.real() : value.real() );
                max_carriers = std::max(max_carriers, carriers[ik * nbnd + ib]);
                min_eig = std::min(min_eig, energies[ik](ib));
                max_eig = std::max(max_eig, energies[ik](ib));
            }
        }

        std::ofstream os(filename__);
        os << "# fat bands of the populations, on the path of kpath (interpolated from the EDUS k grid)\n";
        os << "# k-number    energy(eV)    population P_n(k)    carriers (1-P for filled bands, P for empty bands)\n";
        for (int ik = 0; ik < num_k; ++ik) {
            for (int ib = 0; ib < nbnd; ++ib) {
                os << std::setw(6) << ik;
                os << std::setw(15) << std::setprecision(6) << Convert(energies[ik](ib), AuEnergy, ElectronVolt);
                os << std::setw(22) << std::setprecision(12) << population[ik * nbnd + ib];
                os << std::setw(22) << std::setprecision(12) << carriers[ik * nbnd + ib];
                os << "\n";
            }
        }

        /* gnuplot script: bands as small points, carriers as circles with size proportional to the population */
        auto gnu = std::filesystem::path(filename__).replace_extension(".gnu");
        std::ofstream Gnuplot(gnu);
        min_eig = Convert(min_eig, AuEnergy, ElectronVolt);
        max_eig = Convert(max_eig, AuEnergy, ElectronVolt);
        auto eig_range = max_eig - min_eig;
        Gnuplot << "set xrange [ " << -double(num_k)/10. << " : " << 11./10.*double(num_k) << "]" << std::endl;
        Gnuplot << "set yrange [" << min_eig - eig_range/10. << ": " << max_eig + eig_range/10. << "]" << std::endl;
        for (int ivertex = 0; ivertex < int(vertices.size()); ++ivertex) {
            int kpt = MeshGridPath.find(vertices[ivertex], ivertex);
            Gnuplot << "set arrow from  " << kpt << ", " << min_eig - eig_range/10. << " to " << kpt << ", " << max_eig + eig_range/10. << " nohead" << std::endl;
        }
        Gnuplot << "set ylabel \"energy (eV)\"" << std::endl;
        Gnuplot << "maxsize = 3." << std::endl;
        Gnuplot << "maxcarriers = " << std::max(max_carriers, 1.e-30) << std::endl;
        Gnuplot << "threshold = 1.e-3*maxcarriers" << std::endl;
        Gnuplot << "plot \"" << std::filesystem::path(filename__).filename().string() << "\" u 1:2 w p pt 7 ps 0.2 lc rgb \"grey\" notitle, \\" << std::endl;
        Gnuplot << "     \"\" u 1:($4 > threshold ? $2 : 1/0):(maxsize*$4/maxcarriers) w p pt 6 ps variable lc rgb \"red\" title \"carriers\"" << std::endl;
        Gnuplot << "pause -1" << std::endl;
        output::print("   fat bands on the path -> ", filename__, " (max carriers ", max_carriers, ")");
    }
} // namespace electron
