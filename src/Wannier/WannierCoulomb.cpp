/// Density-density matrix elements of the bare and of the statically screened Coulomb interaction between
/// Wannier functions,
/// @f[ U_{nm}(R) = \int\!\!\int dr\, dr'\, |w_{n0}(r)|^2\, W(r,r')\, |w_{mR}(r')|^2
///              = \frac{\Omega}{N_k^3}\sum_q e^{-iq\cdot R} \sum_{GG'} \hat\rho^{q*}_{nn}(G)\, W_{GG'}(q)\,
///                \hat\rho^{q}_{mm}(G') @f]
/// with @f$ \hat\rho^q_{nn}(G) = \frac{1}{\Omega}\int_\Omega dr\, e^{-iG\cdot r} \sum_k u^*_{nk}(r) u_{n,k+q}(r) @f$
/// (the pair density of WannierWavefunctions at R = 0) and
/// @f$ W_{GG'}(q) = \frac{4\pi}{|q+G|\,|q+G'|}\,\tilde\epsilon^{-1}_{GG'}(q) @f$ (Hartree atomic units).
/// The bare interaction has @f$ \tilde\epsilon^{-1} = \delta_{GG'} @f$. The screening is read from the ndb.em1s
/// database of yambo/Lumen (yambo -X s), where X_Q_<iq> holds the symmetrized
/// @f$ \tilde\epsilon^{-1}_{GG'}(q) - \delta_{GG'} = v^{1/2}(q+G)\,\chi_{GG'}(q)\,v^{1/2}(q+G') @f$
/// (see src/pol_function/X_redux.F of yambo). The q points of the database must cover the whole k grid of
/// wannier90 (yambo run on the nscf with nosym/noinv used for wannier90): no symmetry unfolding is done.
///
/// Output, in the format of the KCW files read by EDUS (ModelCoulomb::initialize_Potential), in Rydberg:
/// <output_dir>/barecoulomb.txt and screencoulomb.txt. For each R of the k grid (from -M/2 to M-1-M/2), a line
/// "R1 R2 R3" followed by the lines "n m ispin Re(U_nm(R)) Im(U_nm(R))", n, m = 1..num_wann, ispin = 1, 2
/// (the same value, as in KCW). n is the Wannier function in the cell 0, m the one in the cell R.
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>

#ifdef EDUS_MPI
#include <fftw3-mpi.h>
#include "core/mpi/Communicator.hpp"
#else
#include <fftw3.h>
#endif

#ifdef EDUS_HDF5
#include <hdf5.h>
#endif

#include "Constants.hpp"
#include "Wannier/WannierWavefunctions.hpp"

namespace {

int mpi_rank()
{
#ifdef EDUS_MPI
    return mpi::Communicator::world().rank();
#else
    return 0;
#endif
}

int mpi_size()
{
#ifdef EDUS_MPI
    return mpi::Communicator::world().size();
#else
    return 1;
#endif
}

/// smallest n >= n_min whose only prime factors are 2, 3, 5, 7
int good_fft_size(int n_min__)
{
    for (int n = std::max(n_min__, 1);; ++n) {
        int m = n;
        for (int p : {2, 3, 5, 7}) {
            while (m % p == 0) {
                m /= p;
            }
        }
        if (m == 1) {
            return n;
        }
    }
}

inline int positive_mod(int a__, int n__)
{
    int r = a__ % n__;
    return r < 0 ? r + n__ : r;
}

#ifdef EDUS_HDF5
/// Reads a whole dataset of a netCDF-4 (HDF5) file of yambo as double. dims__: dimensions in C order
/// (the reverse of the Fortran order of yambo)
std::vector<double> h5_read(const std::string& file__, const std::string& name__, std::vector<hsize_t>& dims__)
{
    hid_t file = H5Fopen(file__.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) {
        throw std::runtime_error("WannierCoulomb: cannot open " + file__);
    }
    if (H5Lexists(file, name__.c_str(), H5P_DEFAULT) <= 0) {
        H5Fclose(file);
        throw std::runtime_error("WannierCoulomb: no variable " + name__ + " in " + file__);
    }
    hid_t dset  = H5Dopen2(file, name__.c_str(), H5P_DEFAULT);
    hid_t space = H5Dget_space(dset);
    int rank    = H5Sget_simple_extent_ndims(space);
    dims__.assign(rank, 0);
    H5Sget_simple_extent_dims(space, dims__.data(), nullptr);
    size_t n = 1;
    for (auto d : dims__) {
        n *= d;
    }
    std::vector<double> out(n);
    herr_t status = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, out.data());
    H5Sclose(space);
    H5Dclose(dset);
    H5Fclose(file);
    if (status < 0) {
        throw std::runtime_error("WannierCoulomb: cannot read " + name__ + " from " + file__);
    }
    return out;
}

bool h5_has(const std::string& file__, const std::string& name__)
{
    if (!std::filesystem::exists(file__)) {
        return false;
    }
    hid_t file = H5Fopen(file__.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) {
        return false;
    }
    bool has = H5Lexists(file, name__.c_str(), H5P_DEFAULT) > 0;
    H5Fclose(file);
    return has;
}
#endif

/// ndb.em1s of yambo/Lumen converted to the crystal coordinates of the QE cell
struct YamboEm1s
{
    std::string dir;
    /// q points (crystal coordinates of the reciprocal lattice) and G vectors (integer, crystal)
    std::vector<std::array<double, 3>> q;
    std::vector<std::array<int, 3>> G;

    /// X_{GG'}(q) = eps~^-1_{GG'}(q) - delta_{GG'} at the first frequency, layout [ig * ng + igp]
    std::vector<std::complex<double>> X(int iq__) const
    {
#ifdef EDUS_HDF5
        const std::string var  = "X_Q_" + std::to_string(iq__ + 1);
        std::string file       = dir + "/ndb.em1s";
        if (!h5_has(file, var)) {
            file = dir + "/ndb.em1s_fragment_" + std::to_string(iq__ + 1);
        }
        std::vector<hsize_t> dims;
        auto data      = h5_read(file, var, dims);
        const size_t n = G.size();
        if (dims.size() != 4 || dims[1] != n || dims[2] != n || dims[3] != 2) {
            throw std::runtime_error("WannierCoulomb: unexpected shape of " + var + " in " + file);
        }
        // C order [iw][j][i][re/im] = Fortran X(re/im, i, j, iw), with i the row (G) and j the column (G')
        std::vector<std::complex<double>> out(n * n);
        for (size_t j = 0; j < n; ++j) {
            for (size_t i = 0; i < n; ++i) {
                out[i * n + j] = {data[(j * n + i) * 2], data[(j * n + i) * 2 + 1]};
            }
        }
        return out;
#else
        (void)iq__;
        return {};
#endif
    }
};

YamboEm1s read_em1s(const std::string& dir__, const std::string& save__, const std::array<std::array<double, 3>, 3>& a__)
{
    YamboEm1s db;
    db.dir = dir__;
#ifdef EDUS_HDF5
    std::vector<hsize_t> dims;
    auto alat = h5_read(save__ + "/ns.db1", "LATTICE_PARAMETER", dims);
    auto qpt  = h5_read(dir__ + "/ndb.em1s", "HEAD_QPT", dims);
    const int nq = int(dims.at(1));
    auto gv   = h5_read(dir__ + "/ndb.em1s", "X_RL_vecs", dims);
    const int ng = int(dims.at(1));
    // yambo units (iku): v_cart_i = v_i 2 pi / alat_i  ->  v_crystal_j = v_cart . a_j / (2 pi) = sum_i v_i a_j[i] / alat_i
    auto to_crystal = [&](const std::vector<double>& v, int n, int i) {
        std::array<double, 3> c{0., 0., 0.};
        for (int j = 0; j < 3; ++j) {
            for (int x = 0; x < 3; ++x) {
                c[j] += v[size_t(x) * n + i] * a__[j][x] / alat[x];
            }
        }
        return c;
    };
    for (int iq = 0; iq < nq; ++iq) {
        db.q.push_back(to_crystal(qpt, nq, iq));
    }
    for (int ig = 0; ig < ng; ++ig) {
        auto c = to_crystal(gv, ng, ig);
        std::array<int, 3> Gi;
        for (int x = 0; x < 3; ++x) {
            Gi[x] = int(std::lround(c[x]));
            if (std::abs(c[x] - Gi[x]) > 1.e-4) {
                std::stringstream ss;
                ss << "WannierCoulomb: G vector " << ig + 1 << " of " << dir__ << "/ndb.em1s has non-integer crystal "
                   << "coordinates (" << c[0] << ", " << c[1] << ", " << c[2] << "): different cell than QE?";
                throw std::runtime_error(ss.str());
            }
        }
        db.G.push_back(Gi);
    }
#else
    (void)save__;
    (void)a__;
    throw std::runtime_error("WannierCoulomb: reading " + dir__ + "/ndb.em1s requires EDUS built with -DEDUS_HDF5=ON");
#endif
    return db;
}

/// ndb.RIM of yambo/Lumen: RIM_qpg(q, G, G') = int_sBZ(q) d^3q'/(2pi)^3 2 / (|q+q'+G| |q+q'+G'|) on the q points and
/// the first ng G vectors of ndb.em1s (both are prefixes of the G list of the yambo SAVE)
struct YamboRim
{
    int nq = 0, ng = 0;
    std::vector<double> data;  // C order [G'][G][q]
    /// int_sBZ(q) d^3q'/(2pi)^3 2/|q+q'+G|^2 (diagonal, yambo default QpgFull = no)
    double diag(int iq__, int ig__) const { return data[(size_t(ig__) * ng + ig__) * nq + iq__]; }
};

YamboRim read_rim(const std::string& file__, const std::string& save__, const std::array<std::array<double, 3>, 3>& a__,
                  const YamboEm1s& db__)
{
    YamboRim r;
#ifdef EDUS_HDF5
    std::vector<hsize_t> dims;
    r.data = h5_read(file__, "RIM_qpg", dims);
    if (dims.size() != 3 || dims[0] != dims[1]) {
        throw std::runtime_error("WannierCoulomb: unexpected shape of RIM_qpg in " + file__);
    }
    r.ng = int(dims[0]);
    r.nq = int(dims[2]);
    if (r.nq != int(db__.q.size()) || r.ng > int(db__.G.size())) {
        throw std::runtime_error("WannierCoulomb: " + file__ + " has a different q grid or more G vectors than ndb.em1s");
    }
    auto alat = h5_read(save__ + "/ns.db1", "LATTICE_PARAMETER", dims);
    auto qpt  = h5_read(file__, "HEAD_QPT", dims);
    for (int iq = 0; iq < r.nq; ++iq) {
        for (int j = 0; j < 3; ++j) {
            double c = 0.;
            for (int x = 0; x < 3; ++x) {
                c += qpt[size_t(x) * r.nq + iq] * a__[j][x] / alat[x];
            }
            if (std::abs(c - db__.q[iq][j]) > 1.e-5) {
                throw std::runtime_error("WannierCoulomb: the q points of " + file__ + " differ from those of ndb.em1s");
            }
        }
    }
#else
    (void)file__;
    (void)save__;
    (void)a__;
    (void)db__;
    throw std::runtime_error("WannierCoulomb: rim_file needs EDUS_HDF5");
#endif
    return r;
}

/// average of |q|^-p (p = 1, 2) on the small cell of the q grid around Gamma, by midpoint sampling
/// ndb.cutoff of yambo/Lumen: CUT_BARE_QPG(G, q) = the (cut) |q+G| that yambo uses everywhere, v(q+G) = 4pi/|q+G|^2,
/// on the G list of the SAVE (ns.db1 "G-VECTORS", same order) and the q points of ndb.em1s. With the slab cutoff and
/// the RIM yambo leaves at q = G = 0 the value regularized at |q| = q0_def_norm. G maps a G vector (integer, crystal)
/// to its index in the list
struct YamboCutoff
{
    int nq = 0, ng = 0;
    std::vector<double> data;  // C order [G][q][re/im]
    std::map<std::array<int, 3>, int> G;
    double qpg(int iq__, int ig__) const
    {
        return std::abs(std::complex<double>(data[(size_t(ig__) * nq + iq__) * 2], data[(size_t(ig__) * nq + iq__) * 2 + 1]));
    }
};

YamboCutoff read_cutoff(const std::string& file__, const std::string& save__, const std::array<std::array<double, 3>, 3>& a__,
                        const YamboEm1s& db__)
{
    YamboCutoff c;
#ifdef EDUS_HDF5
    std::vector<hsize_t> dims;
    c.data = h5_read(file__, "CUT_BARE_QPG", dims);
    if (dims.size() != 3 || dims[2] != 2) {
        throw std::runtime_error("WannierCoulomb: unexpected shape of CUT_BARE_QPG in " + file__);
    }
    c.ng = int(dims[0]);
    c.nq = int(dims[1]);
    if (c.nq != int(db__.q.size())) {
        throw std::runtime_error("WannierCoulomb: " + file__ + " has a different q grid than ndb.em1s");
    }
    auto alat = h5_read(save__ + "/ns.db1", "LATTICE_PARAMETER", dims);
    auto qpt  = h5_read(file__, "HEAD_QPT", dims);
    for (int iq = 0; iq < c.nq; ++iq) {
        for (int j = 0; j < 3; ++j) {
            double x = 0.;
            for (int i = 0; i < 3; ++i) {
                x += qpt[size_t(i) * c.nq + iq] * a__[j][i] / alat[i];
            }
            if (std::abs(x - db__.q[iq][j]) > 1.e-5) {
                throw std::runtime_error("WannierCoulomb: the q points of " + file__ + " differ from those of ndb.em1s");
            }
        }
    }
    auto gv = h5_read(save__ + "/ns.db1", "G-VECTORS", dims);
    const int ngall = int(dims.at(1));
    if (ngall < c.ng) {
        throw std::runtime_error("WannierCoulomb: " + file__ + " has more G vectors than the SAVE");
    }
    for (int ig = 0; ig < c.ng; ++ig) {
        std::array<int, 3> Gi;
        for (int j = 0; j < 3; ++j) {
            double x = 0.;
            for (int i = 0; i < 3; ++i) {
                x += gv[size_t(i) * ngall + ig] * a__[j][i] / alat[i];
            }
            Gi[j] = int(std::lround(x));
        }
        c.G[Gi] = ig;
    }
#else
    (void)file__;
    (void)save__;
    (void)a__;
    (void)db__;
    throw std::runtime_error("WannierCoulomb: reading ndb.cutoff needs EDUS_HDF5");
#endif
    return c;
}

/// average of f(q) on the small cell of the q grid around Gamma (q_j in [-1/2, 1/2] b_j / M_j); with slab__ only
/// the in-plane cell (q_3 = 0)
template <typename F>
double mini_bz_average(const std::array<std::array<double, 3>, 3>& b__, const std::array<int, 3>& M__, const bool slab__,
                       F&& f__)
{
    const int Ns  = 96;
    const int Ns3 = slab__ ? 1 : Ns;
    double s      = 0.;
    for (int i1 = 0; i1 < Ns; ++i1) {
        for (int i2 = 0; i2 < Ns; ++i2) {
            for (int i3 = 0; i3 < Ns3; ++i3) {
                double x[3] = {(i1 + 0.5) / Ns - 0.5, (i2 + 0.5) / Ns - 0.5, slab__ ? 0. : (i3 + 0.5) / Ns - 0.5};
                std::array<double, 3> q{0., 0., 0.};
                for (int c = 0; c < 3; ++c) {
                    for (int j = 0; j < 3; ++j) {
                        q[c] += x[j] * b__[j][c] / M__[j];
                    }
                }
                s += f__(q);
            }
        }
    }
    return s / (double(Ns) * Ns * Ns3);
}

} // namespace

void WannierWavefunctions::compute_coulomb()
{
    const auto& cp   = p_.coulomb;
    const int nk     = gauge_.num_kpts();
    const int nw     = gauge_.num_wann();
    const auto& M    = gauge_.mp_grid();
    const int rank   = mpi_rank();
    const bool screen = !cp.em1s_dir.empty();

    auto cart = [&](const std::array<double, 3>& v) {
        std::array<double, 3> c{0., 0., 0.};
        for (int j = 0; j < 3; ++j) {
            for (int x = 0; x < 3; ++x) {
                c[x] += v[j] * b_[j][x];
            }
        }
        return c;
    };
    auto norm2 = [](const std::array<double, 3>& c) { return c[0] * c[0] + c[1] * c[1] + c[2] * c[2]; };
    // Coulomb kernel v(Q) (Hartree a.u.), with the slab cutoff of yambo for 2D systems
    const bool slab = (cp.cutoff == "slab");
    const double Lslab = std::sqrt(norm2(a_[2]));
    if (slab && (std::abs(a_[2][0]) + std::abs(a_[2][1]) > 1.e-8 || std::abs(a_[0][2]) + std::abs(a_[1][2]) > 1.e-8)) {
        throw std::runtime_error("WannierCoulomb: the slab cutoff needs a3 along z and a1, a2 in the plane");
    }
    auto vcoul = [&](const std::array<double, 3>& Q) {
        double Q2 = norm2(Q);
        double f  = 1.;
        if (slab) {
            f = 1. - std::exp(-std::sqrt(Q[0] * Q[0] + Q[1] * Q[1]) * Lslab / 2.) * std::cos(Q[2] * Lslab / 2.);
        }
        return 4. * pi * f / Q2;
    };

    // ---- q points and G vectors -------------------------------------------------------------------------------
    YamboEm1s db;
    std::vector<std::array<double, 3>> qlist;
    if (screen) {
        db    = read_em1s(cp.em1s_dir, cp.yambo_save, a_);
        qlist = db.q;
    } else {
        if (cp.ecut_bare <= 0.) {
            throw std::runtime_error("WannierCoulomb: without em1s_dir, ecut_bare (Ry) must be given");
        }
        for (int i1 = 0; i1 < M[0]; ++i1) {
            for (int i2 = 0; i2 < M[1]; ++i2) {
                for (int i3 = 0; i3 < M[2]; ++i3) {
                    qlist.push_back({double(i1) / M[0], double(i2) / M[1], double(i3) / M[2]});
                }
            }
        }
    }
    // yambo RIM: "head" -> only the q -> 0 head and wings; "full" -> as yambo: every v(q+G) with G among the first
    // RIM_ng G vectors replaced by its average on the cell of q (diagonal), the G = 0 wings by <v(q)>^1/2
    const bool rim_full = !cp.rim_file.empty() && cp.rim_mode == "full";
    YamboRim rim;
    if (!cp.rim_file.empty()) {
        if (!screen || (!slab && !cp.include_head)) {
            throw std::runtime_error("WannierCoulomb: rim_file needs em1s_dir, and include_head without cutoff");
        }
        if (cp.rim_mode != "head" && cp.rim_mode != "full") {
            throw std::runtime_error("WannierCoulomb: rim_mode must be head or full");
        }
        if (slab && cp.rim_mode != "full") {
            throw std::runtime_error("WannierCoulomb: with cutoff = slab the RIM is used as yambo does, rim_mode = full");
        }
        rim = read_rim(cp.rim_file, cp.yambo_save, a_, db);
    }
    // slab + RIM: exactly scatter_Gamp of yambo. v(q+G) = 4pi/|q+G|^2 with the cut |q+G| of yambo (ndb.cutoff, next to
    // ndb.em1s) for every G, q = G = 0 included (regularized by yambo); the RIM average on the diagonal of the RIM G
    // vectors (head included); no other treatment of head and wings (yambo skips it with a cutoff)
    const bool yambo_slab = slab && !cp.rim_file.empty();
    YamboCutoff ycut;
    if (yambo_slab) {
        ycut = read_cutoff(cp.em1s_dir + "/ndb.cutoff", cp.yambo_save, a_, db);
    }
    double max_vcut_diff = 0.;
    const int nq = int(qlist.size());
    // every point of the k grid exactly once (modulo G)
    {
        std::vector<int> seen(size_t(M[0]) * M[1] * M[2], 0);
        for (auto& q : qlist) {
            int idx[3];
            for (int x = 0; x < 3; ++x) {
                double v = q[x] * M[x];
                idx[x]   = int(std::lround(v));
                if (std::abs(v - idx[x]) > 1.e-4) {
                    std::stringstream ss;
                    ss << "WannierCoulomb: q = (" << q[0] << ", " << q[1] << ", " << q[2]
                       << ") of ndb.em1s is not on the wannier90 k grid (mp_grid " << M[0] << " " << M[1] << " " << M[2]
                       << ")";
                    throw std::runtime_error(ss.str());
                }
                idx[x] = positive_mod(idx[x], M[x]);
            }
            seen[(size_t(idx[0]) * M[1] + idx[1]) * M[2] + idx[2]]++;
        }
        for (auto s : seen) {
            if (s != 1) {
                throw std::runtime_error(
                    "WannierCoulomb: the q points of ndb.em1s do not cover the k grid exactly once. Run yambo on the "
                    "nscf used for wannier90 (nosym, noinv), so that no symmetry unfolding is needed");
            }
        }
    }
    // for each q: block of the screening (indices in db.G) and G vectors of the bare interaction
    // Gs[iq]: all the G vectors used; nblock: the first nblock are those of ndb.em1s, in the same order
    const int nblock = screen ? int(db.G.size()) : 0;
    std::vector<std::vector<std::array<int, 3>>> Gs(nq);
    int Gmax[3] = {0, 0, 0};
    for (int iq = 0; iq < nq; ++iq) {
        Gs[iq] = db.G;
        if (cp.ecut_bare > 0.) {
            // sphere |q+G|^2 <= ecut_bare (Ry, bohr^-2), on top of the G vectors of the screening block
            double bmin = 1.e30;
            for (int j = 0; j < 3; ++j) {
                bmin = std::min(bmin, std::sqrt(norm2(b_[j])));
            }
            // |G_j| <= |q+G| |a_j| / (2 pi) + |q_j|
            int nmax[3];
            for (int j = 0; j < 3; ++j) {
                nmax[j] = int(std::sqrt(cp.ecut_bare) * std::sqrt(norm2(a_[j])) / (2. * pi) + std::abs(qlist[iq][j])) + 1;
            }
            std::map<std::array<int, 3>, int> in_block;
            for (auto& G : db.G) {
                in_block[G] = 1;
            }
            for (int g1 = -nmax[0]; g1 <= nmax[0]; ++g1) {
                for (int g2 = -nmax[1]; g2 <= nmax[1]; ++g2) {
                    for (int g3 = -nmax[2]; g3 <= nmax[2]; ++g3) {
                        std::array<int, 3> G{g1, g2, g3};
                        auto Q = cart({qlist[iq][0] + g1, qlist[iq][1] + g2, qlist[iq][2] + g3});
                        if (norm2(Q) <= cp.ecut_bare && !in_block.count(G)) {
                            Gs[iq].push_back(G);
                        }
                    }
                }
            }
        }
        for (auto& G : Gs[iq]) {
            for (int x = 0; x < 3; ++x) {
                Gmax[x] = std::max(Gmax[x], std::abs(G[x]));
            }
        }
    }

    // ---- k + q = k_j + Gshift ---------------------------------------------------------------------------------
    std::vector<std::vector<int>> ikq(nq, std::vector<int>(nk, -1));
    std::vector<std::vector<std::array<int, 3>>> Gshift(nq, std::vector<std::array<int, 3>>(nk));
    int Smax[3] = {0, 0, 0};
    for (int iq = 0; iq < nq; ++iq) {
        for (int ik = 0; ik < nk; ++ik) {
            auto k = gauge_.kpt(ik);
            for (int jk = 0; jk < nk && ikq[iq][ik] < 0; ++jk) {
                auto kj = gauge_.kpt(jk);
                std::array<int, 3> G;
                bool match = true;
                for (int x = 0; x < 3; ++x) {
                    double d = k[x] + qlist[iq][x] - kj[x];
                    G[x]     = int(std::lround(d));
                    match    = match && std::abs(d - G[x]) < 1.e-4;
                }
                if (match) {
                    ikq[iq][ik]    = jk;
                    Gshift[iq][ik] = G;
                    for (int x = 0; x < 3; ++x) {
                        Smax[x] = std::max(Smax[x], std::abs(G[x]));
                    }
                }
            }
            if (ikq[iq][ik] < 0) {
                throw std::runtime_error("WannierCoulomb: k + q not found in the wannier90 k list");
            }
        }
    }

    // ---- Bloch functions in the Wannier gauge and FFT grid ------------------------------------------------------
    std::vector<QEWavefunction> wk(nk);
    int max_miller[3] = {0, 0, 0};
    for (int ik = 0; ik < nk; ++ik) {
        wk[ik] = bloch_wannier_gauge(ik);
        for (int ig = 0; ig < wk[ik].igwx; ++ig) {
            for (int x = 0; x < 3; ++x) {
                max_miller[x] = std::max(max_miller[x], std::abs(wk[ik].miller(ig, x)));
            }
        }
    }
    // u*_k u_{k_j} e^{-iGshift.r} has |m| <= 2 max_miller + |Gshift|: no aliasing on |G| <= Gmax if N > that + Gmax
    std::array<int, 3> grid;
    for (int x = 0; x < 3; ++x) {
        grid[x] = good_fft_size(std::max(2 * max_miller[x] + Smax[x] + Gmax[x] + 1, 2 * Gmax[x] + 1));
    }
    const size_t Nr = size_t(grid[0]) * grid[1] * grid[2];
    std::vector<int> bands(nw);
    for (int n = 0; n < nw; ++n) {
        bands[n] = n;
    }

    if (rank == 0) {
        std::cout << "\n  Coulomb matrix elements U_nm(R) (density-density, Wannier basis)\n"
                  << "    screening        : " << (screen ? cp.em1s_dir + "/ndb.em1s" : std::string("none (bare only)"))
                  << "\n    cutoff           : " << cp.cutoff
                  << "\n    q points         : " << nq << "\n    G vectors        : " << nblock << " (screening block), "
                  << Gs[0].size() << " (bare, q = 1)\n    FFT grid (cell)  : " << grid[0] << " " << grid[1] << " "
                  << grid[2] << "\n    memory per rank  : " << std::fixed << std::setprecision(1)
                  << double(nk) * nw * npol_ * Nr * sizeof(std::complex<double>) / double(1 << 20) << " MB (u_nk(r))"
                  << std::defaultfloat << "\n";
    }

    // u_nk(r) for all k, on every rank
    std::vector<std::vector<std::complex<double>>> uk(nk);
    for (int ik = 0; ik < nk; ++ik) {
        periodic_part(wk[ik], {0, 0, 0}, grid, bands, uk[ik]);
    }
    wk.clear();

    // ---- per q: rho^q_nn(G), V_q and W_q --------------------------------------------------------------------------
    // Vq[iq][n][m] = sum_G rho*_n v rho_m (the q + G = 0 term excluded, added later if include_head)
    std::vector<std::complex<double>> Vq(size_t(nq) * nw * nw, 0.), Wq(size_t(nq) * nw * nw, 0.);
    // head and wings at q = 0: rho_n(G=0) and the wings sum_G' X_{0G'} / |G'| rho_m(G') (and its mirror)
    int iq0 = -1;
    std::vector<std::complex<double>> rho0(nw, 0.), wingL(nw, 0.), wingR(nw, 0.);
    std::complex<double> X00 = 0.;
    double max_herm = 0.;

    auto* buffer = reinterpret_cast<std::complex<double>*>(fftw_malloc(sizeof(fftw_complex) * Nr));
    auto plan    = fftw_plan_dft_3d(grid[0], grid[1], grid[2], reinterpret_cast<fftw_complex*>(buffer),
                                    reinterpret_cast<fftw_complex*>(buffer), FFTW_FORWARD, FFTW_ESTIMATE);
    std::vector<std::complex<double>> ph1(grid[0]), ph2(grid[1]), ph3(grid[2]);
    for (int iq = rank; iq < nq; iq += mpi_size()) {
        const auto& q  = qlist[iq];
        const bool isq0 = std::abs(q[0] - std::round(q[0])) + std::abs(q[1] - std::round(q[1])) +
                              std::abs(q[2] - std::round(q[2])) < 1.e-8;
        const int ng   = int(Gs[iq].size());
        // rho_n(G) for the G of this q, layout [n * ng + ig]
        std::vector<std::complex<double>> rho(size_t(nw) * ng);
        for (int n = 0; n < nw; ++n) {
            std::fill(buffer, buffer + Nr, std::complex<double>(0.));
            for (int ik = 0; ik < nk; ++ik) {
                const int jk  = ikq[iq][ik];
                const auto& G = Gshift[iq][ik];
                // u_{k+q}(r) = e^{-iG.r} u_{k_j}(r)
                for (int i = 0; i < grid[0]; ++i) {
                    ph1[i] = std::exp(-im * 2. * pi * double(G[0] * i) / double(grid[0]));
                }
                for (int i = 0; i < grid[1]; ++i) {
                    ph2[i] = std::exp(-im * 2. * pi * double(G[1] * i) / double(grid[1]));
                }
                for (int i = 0; i < grid[2]; ++i) {
                    ph3[i] = std::exp(-im * 2. * pi * double(G[2] * i) / double(grid[2]));
                }
                for (int ipol = 0; ipol < npol_; ++ipol) {
                    const auto* ua = uk[ik].data() + (size_t(n) * npol_ + ipol) * Nr;
                    const auto* ub = uk[jk].data() + (size_t(n) * npol_ + ipol) * Nr;
                    #pragma omp parallel for collapse(2)
                    for (int i1 = 0; i1 < grid[0]; ++i1) {
                        for (int i2 = 0; i2 < grid[1]; ++i2) {
                            const auto p12 = ph1[i1] * ph2[i2];
                            const size_t o = (size_t(i1) * grid[1] + i2) * grid[2];
                            for (int i3 = 0; i3 < grid[2]; ++i3) {
                                buffer[o + i3] += std::conj(ua[o + i3]) * ub[o + i3] * p12 * ph3[i3];
                            }
                        }
                    }
                }
            }
            fftw_execute(plan);  // sum_r rho(r) e^{-iG.r}
            for (int ig = 0; ig < ng; ++ig) {
                const auto& G = Gs[iq][ig];
                size_t idx    = (size_t(positive_mod(G[0], grid[0])) * grid[1] + positive_mod(G[1], grid[1])) * grid[2] +
                             positive_mod(G[2], grid[2]);
                rho[size_t(n) * ng + ig] = buffer[idx] / double(Nr);
            }
        }

        // |q+G|, and the index of q + G = 0
        std::vector<double> qG(ng), vG(ng, 0.), svG(ng, 0.);
        int ig0 = -1;
        for (int ig = 0; ig < ng; ++ig) {
            const auto& G = Gs[iq][ig];
            auto Q        = cart({q[0] + G[0], q[1] + G[1], q[2] + G[2]});
            qG[ig]        = std::sqrt(norm2(Q));
            if (qG[ig] < 1.e-8) {
                ig0 = ig;
            } else {
                vG[ig]  = vcoul(Q);
                svG[ig] = std::sqrt(vG[ig]);
            }
            if (yambo_slab) {
                auto it = ycut.G.find(G);
                if (it == ycut.G.end()) {
                    throw std::runtime_error("WannierCoulomb: a G vector of the bare interaction is not in ndb.cutoff");
                }
                const double qpg = ycut.qpg(iq, it->second);
                const double vy  = 4. * pi / (qpg * qpg);
                // (the slab v vanishes for q + G along z with cos(Q_z L/2) = 1, where yambo regularizes it)
                if (qG[ig] >= 1.e-8 && vG[ig] > 1.e-8 * vy) {
                    max_vcut_diff = std::max(max_vcut_diff, std::abs(vy - vG[ig]) / vG[ig]);
                }
                vG[ig]  = vy;
                svG[ig] = std::sqrt(vy);
            }
        }
        // with yambo_slab the term q + G = 0 is a G vector like the others
        if (yambo_slab) {
            ig0 = -1;
        }
        // RIM: vG = <v(q+G)> on the cell of q (diagonal terms), svG(G=0) = <v(q)>^1/2 (wings, q != 0); the
        // off-diagonal G != G' keep the point values v^1/2(q+G) v^1/2(q+G'), as scatter_Gamp of yambo
        std::vector<double> vD = vG;
        if (rim_full) {
            for (int ig = 0; ig < rim.ng; ++ig) {
                if (ig == ig0) {
                    continue;
                }
                vD[ig] = 4. * pi * rim.diag(iq, ig) * omega_ * nk / 2.;
                const auto& G = Gs[iq][ig];
                if (!yambo_slab && G[0] == 0 && G[1] == 0 && G[2] == 0) {
                    svG[ig] = std::sqrt(vD[ig]);
                }
            }
        }
        std::vector<std::complex<double>> X;
        if (screen) {
            X = db.X(iq);
            for (int i = 0; i < nblock; ++i) {
                for (int j = 0; j < nblock; ++j) {
                    max_herm = std::max(max_herm, std::abs(X[i * nblock + j] - std::conj(X[j * nblock + i])));
                }
            }
        }
        for (int n = 0; n < nw; ++n) {
            for (int m = 0; m < nw; ++m) {
                const auto* rn = rho.data() + size_t(n) * ng;
                const auto* rm = rho.data() + size_t(m) * ng;
                std::complex<double> v = 0., w = 0.;
                for (int ig = 0; ig < ng; ++ig) {
                    if (ig == ig0) {
                        continue;
                    }
                    v += std::conj(rn[ig]) * rm[ig] * vD[ig];
                }
                w = v;
                if (screen) {
                    // + sum_{GG' in block} rho*_n(G) v^1/2(q+G) X_{GG'} v^1/2(q+G') rho_m(G')
                    for (int i = 0; i < nblock; ++i) {
                        if (i == ig0) {
                            continue;
                        }
                        std::complex<double> s = 0.;
                        for (int j = 0; j < nblock; ++j) {
                            if (j == ig0) {
                                continue;
                            }
                            s += X[size_t(i) * nblock + j] * rm[j] * svG[j];
                        }
                        w += std::conj(rn[i]) * svG[i] * s;
                        if (rim_full) {
                            // diagonal G = G': <v(q+G)> instead of v^1/2 v^1/2
                            w += std::conj(rn[i]) * X[size_t(i) * nblock + i] * rm[i] * (vD[i] - svG[i] * svG[i]);
                        }
                    }
                }
                Vq[(size_t(iq) * nw + n) * nw + m] = v;
                Wq[(size_t(iq) * nw + n) * nw + m] = w;
            }
        }
        if (isq0 && !yambo_slab) {
            if (ig0 < 0) {
                throw std::runtime_error("WannierCoulomb: G = -q missing among the G vectors at q = 0");
            }
            iq0 = iq;
            for (int n = 0; n < nw; ++n) {
                rho0[n] = rho[size_t(n) * ng + ig0];
                if (screen && ig0 < nblock) {
                    for (int j = 0; j < nblock; ++j) {
                        if (j != ig0) {
                            wingL[n] += X[size_t(ig0) * nblock + j] * rho[size_t(n) * ng + j] * svG[j];
                            wingR[n] += std::conj(rho[size_t(n) * ng + j]) * X[size_t(j) * nblock + ig0] * svG[j];
                        }
                    }
                }
            }
            if (screen && ig0 < nblock) {
                X00 = X[size_t(ig0) * nblock + ig0];
            }
        }
    }
    fftw_destroy_plan(plan);
    fftw_free(buffer);

#ifdef EDUS_MPI
    {
        auto comm = mpi::Communicator::world().communicator();
        auto sum  = [&](std::vector<std::complex<double>>& v) {
            MPI_Allreduce(MPI_IN_PLACE, v.data(), int(v.size()), MPI_C_DOUBLE_COMPLEX, MPI_SUM, comm);
        };
        sum(Vq);
        sum(Wq);
        sum(rho0);
        sum(wingL);
        sum(wingR);
        std::vector<std::complex<double>> x00{X00};
        sum(x00);
        X00 = x00[0];
        MPI_Allreduce(MPI_IN_PLACE, &iq0, 1, MPI_INT, MPI_MAX, comm);
        MPI_Allreduce(MPI_IN_PLACE, &max_herm, 1, MPI_DOUBLE, MPI_MAX, comm);
        MPI_Allreduce(MPI_IN_PLACE, &max_vcut_diff, 1, MPI_DOUBLE, MPI_MAX, comm);
    }
#endif
    if (rank != 0) {
        return;
    }

    // ---- head and wings (q -> 0), averaged on the small cell around Gamma -------------------------------------------
    if (yambo_slab) {
        std::cout << std::setprecision(6) << "    yambo slab + RIM : v(q+G) from " << cp.em1s_dir << "/ndb.cutoff, RIM ("
                  << rim.ng << " G) on the diagonal, no separate head (as scatter_Gamp); max |v_yambo - v_slab|/v_slab"
                  << " for q + G != 0 = " << max_vcut_diff << "\n";
    }
    if (cp.include_head && !yambo_slab) {
        // <v(q)> and <v^1/2(q)> on the small cell
        double avg2 = mini_bz_average(b_, M, slab, [&](const std::array<double, 3>& Q) { return vcoul(Q); });
        double avg1 = mini_bz_average(b_, M, slab, [&](const std::array<double, 3>& Q) { return std::sqrt(vcoul(Q)); });
        if (!cp.rim_file.empty()) {
            // head <v(q)> = 4pi <1/q^2> from the RIM; wings as yambo (scatter_Gamp): <1/q> ~ R_sphere/2 <1/q^2>,
            // R_sphere the radius of the sphere with the volume of the small cell
            const double inv_q2  = rim.diag(iq0, 0) * omega_ * nk / 2.;
            const double vcell   = std::pow(2. * pi, 3) / (omega_ * nk);
            const double Rsphere = std::cbrt(3. * vcell / (4. * pi));
            std::cout << std::setprecision(6) << "    yambo RIM        : " << cp.rim_file << " (" << cp.rim_mode << ", "
                      << rim.ng << " G), <1/q^2> = " << inv_q2 << " (small cell: " << avg2 / (4. * pi)
                      << "), <1/q> = R_sphere/2 <1/q^2> = " << Rsphere / 2. * inv_q2
                      << " (small cell: " << avg1 / std::sqrt(4. * pi) << ")\n";
            avg2 = 4. * pi * inv_q2;
            avg1 = std::sqrt(4. * pi) * Rsphere / 2. * inv_q2;
        }
        for (int n = 0; n < nw; ++n) {
            for (int m = 0; m < nw; ++m) {
                auto head = std::conj(rho0[n]) * rho0[m] * avg2;
                Vq[(size_t(iq0) * nw + n) * nw + m] += head;
                Wq[(size_t(iq0) * nw + n) * nw + m] +=
                    head * (1. + X00) + avg1 * (std::conj(rho0[n]) * wingL[m] + wingR[n] * rho0[m]);
            }
        }
        std::cout << std::setprecision(6) << "    head included    : <v(q)> = " << avg2 << ", <v^1/2(q)> = " << avg1
                  << " on the small cell, eps~^-1_00(q->0) = " << (1. + X00).real() << "\n";
    }

    // ---- U_nm(R) = Omega / N_k^3 sum_q e^{-iq.R} V_q(n,m), written in Rydberg ------------------------------------------
    const double pref = omega_ / (double(nk) * nk * nk);
    double max_imag   = 0.;
    auto write = [&](const std::vector<std::complex<double>>& Uq, const std::string& name, std::vector<double>& onsite) {
        std::ofstream f(p_.output_dir + "/" + name);
        f << std::scientific << std::setprecision(16);
        onsite.assign(nw, 0.);
        for (int r1 = -M[0] / 2; r1 < M[0] - M[0] / 2; ++r1) {
            for (int r2 = -M[1] / 2; r2 < M[1] - M[1] / 2; ++r2) {
                for (int r3 = -M[2] / 2; r3 < M[2] - M[2] / 2; ++r3) {
                    f << std::setw(7) << r1 << std::setw(7) << r2 << std::setw(7) << r3 << "\n";
                    std::vector<std::complex<double>> U(size_t(nw) * nw, 0.);
                    for (int iq = 0; iq < nq; ++iq) {
                        const auto& q = qlist[iq];
                        auto ph = std::exp(-im * 2. * pi * (q[0] * r1 + q[1] * r2 + q[2] * r3));
                        for (size_t i = 0; i < U.size(); ++i) {
                            U[i] += ph * Uq[size_t(iq) * nw * nw + i];
                        }
                    }
                    for (int n = 0; n < nw; ++n) {
                        for (int m = 0; m < nw; ++m) {
                            auto u = pref * U[size_t(n) * nw + m] * 2.;  // Hartree -> Rydberg
                            max_imag = std::max(max_imag, std::abs(u.imag()));
                            if (r1 == 0 && r2 == 0 && r3 == 0 && n == m) {
                                onsite[n] = u.real();
                            }
                            for (int ispin = 1; ispin <= 2; ++ispin) {
                                f << std::setw(5) << n + 1 << std::setw(5) << m + 1 << std::setw(5) << ispin
                                  << std::setw(28) << u.real() << std::setw(28) << u.imag() << "\n";
                            }
                        }
                    }
                }
            }
        }
    };
    std::vector<double> onV, onW;
    write(Vq, "barecoulomb.txt", onV);
    if (screen) {
        write(Wq, "screencoulomb.txt", onW);
    }
    std::cout << "    on-site U_nn(0) (eV), bare" << (screen ? " / screened" : "") << ":\n";
    for (int n = 0; n < nw; ++n) {
        std::cout << "      " << std::setw(4) << n + 1 << std::fixed << std::setprecision(4) << std::setw(12)
                  << onV[n] * Rydberg_value / ElectronVolt_value;
        if (screen) {
            std::cout << std::setw(12) << onW[n] * Rydberg_value / ElectronVolt_value;
        }
        std::cout << std::defaultfloat << "\n";
    }
    std::cout << "    max |Im U| = " << std::scientific << std::setprecision(3) << max_imag << " Ry";
    if (screen) {
        std::cout << ", max |X - X^+| = " << max_herm;
    }
    std::cout << std::defaultfloat << (cp.include_head ? "" : "  (q + G = 0 term left out, as in KCW)") << "\n"
              << "    written to " << p_.output_dir << "/barecoulomb.txt" << (screen ? " and screencoulomb.txt" : "")
              << " (Rydberg, KCW format read by EDUS)\n";
}
