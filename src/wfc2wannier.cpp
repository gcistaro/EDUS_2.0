/// EDUS_wfc2wannier: reads the Kohn-Sham wavefunctions of Quantum ESPRESSO and rotates them
/// to the Wannier gauge using the matrices written by wannier90 (write_u_matrices = .true.).
///
/// Usage: EDUS_wfc2wannier input.json
///
/// {
///   "qe_save_dir"  : "out/si.save",      // required, contains wfc<ik>.dat
///   "seedname"     : "wannier/si",       // required, <seedname>_u.mat [, _u_dis.mat, .eig]
///   "exclude_bands": [1,2,3,4,5],        // as in the .win file (1-based)
///   "dis_win_min"  : -100.0,             // outer window (eV), only if used in the .win file
///   "dis_win_max"  : 20.0,
///   "spin"         : "none",             // none | up | down (LSDA)
///   "fft_grid"     : [0,0,0],            // 0 = automatic
///   "supercell"    : [3,3,3],
///   "wannier_list" : [],                 // 1-based, empty = all
///   "write_bloch"  : true,               // <output_dir>/wfc<ik>.dat in QE format, num_wann bands
///   "real_space"   : true,               // build w_n(r) on the supercell
///   "write_xsf"    : true,
///   "fix_phase"    : true,
///   "output_dir"   : "wannier_wfc",
///   "pair_density" : {                   // optional: rho^{q,R}_{nn'}(r) = sum_k e^{ik.R} u*_{nk}(r) u_{n'k+q}(r)
///       "q_list"      : [[0,0,0]],       // crystal coordinates, on the k grid
///       "R_list"      : [[0,0,0]],       // crystal coordinates
///       "wannier_list": [],              // 1-based, empty = all
///       "fft_grid"    : [0,0,0]          // 0 = automatic (no aliasing of the product)
///   },                                   // -> <output_dir>/pair_density_q<iq>.dat
///   "coulomb"      : {                   // optional: density-density U_nm(R) of the bare/screened interaction
///       "em1s_dir"    : "em1s",          // yambo/Lumen ndb.em1s (yambo -X s); "" -> bare only. Needs EDUS_HDF5
///       "yambo_save"  : "SAVE",          // yambo SAVE (ns.db1)
///       "ecut_bare"   : 0.,              // Ry, G sphere of the bare interaction; 0 -> G vectors of ndb.em1s
///       "include_head": false,           // add the q -> 0 head and wings (default: left out, as in KCW)
///       "rim_file"    : "",              // yambo ndb.RIM (with include_head): q -> 0 terms as yambo
///       "rim_mode"    : "full"           // "full": RIM on all q and the RIM G vectors (as yambo); "head": only q -> 0
///   }                                    // -> <output_dir>/barecoulomb.txt, screencoulomb.txt (KCW format)
/// }
#include <fstream>
#include <iostream>

#include "Json/json.hpp"
#include "Wannier/WannierWavefunctions.hpp"
#include "core/print_timing.hpp"
#include "core/profiler.hpp"
#include "initialize.hpp"

int main(int argc, char* argv[])
{
    if (argc != 2) {
        std::cerr << "Usage: EDUS_wfc2wannier <input.json>\n";
        return 1;
    }
    initialize();
    {
        PROFILE("wfc2wannier");
        std::ifstream f(argv[1]);
        if (!f.is_open()) {
            throw std::runtime_error(std::string("cannot open ") + argv[1]);
        }
        auto json   = nlohmann::json::parse(f);
        auto params = WannierWavefunctionsParameters::from_json(json);
        WannierWavefunctions wwf(params);
        wwf.run();
    }
    finalize();
    return 0;
}
