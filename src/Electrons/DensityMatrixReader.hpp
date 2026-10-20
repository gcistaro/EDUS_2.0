#ifndef DENSITY_MATRIX_READER_HPP
#define DENSITY_MATRIX_READER_HPP

#include <string>
#include <vector>
#include "Operator/Operator.hpp"
#include "GridStructure/GridStructure.hpp"
#include "Electrons/BandStructure.hpp"
#include "core/mpi/Communicator.hpp"

namespace electron
{
    /// @brief Reads density matrices written by another code, in the wannier gauge, on a k grid.
    /// The directory contains:
    /// - P_in_W_{i1}{i2}_{it}.txt, for i1 <= i2 (0-based, no separator): one line "re im" per k point, with
    ///   @f$ \rho_{i2,i1}(\textbf{k}) @f$ at the snapshot it (the other triangle is fixed by hermiticity)
    /// - the k points of each line, in crystal coordinates (any order and any interval, e.g. [-0.5,0.5))
    /// - optionally the time of each snapshot, one per line
    /// The k points must be the ones of the EDUS grid: each one is mapped to the EDUS index.
    /// Rank 0 reads, the other ranks get their local k points.
    class DensityMatrixReader
    {
        private:
            std::string directory_;
            int num_wann_ = 0;
            /// EDUS global k index of each line of the files
            std::vector<int> edus_index_;
            /// File name (without _it.txt) for each pair i1 <= i2
            std::vector<std::vector<std::string>> basename_;
            int num_snapshots_ = 0;
            /// Time of each snapshot, as written in the file (empty if not found)
            std::vector<double> times_;
            std::string file(int i1__, int i2__, int isnap__) const;

        public:
            DensityMatrixReader() = default;
            /// Reads the k points, maps them to the EDUS grid and counts the snapshots
            void initialize(const std::string& directory__, const std::string& kpoints_file__,
                            const std::string& times_file__, int num_wann__, const MeshGrid& kgrid__,
                            const mpi::Communicator& comm__);
            /// Reads the snapshot isnap__ in DM__ (k space, wannier gauge, local k points)
            void read(int isnap__, Operator<std::complex<double>>& DM__, const MPIindex<3>& mpindex__,
                      const mpi::Communicator& comm__) const;

            int num_snapshots() const { return num_snapshots_; }
            const std::vector<double>& times() const { return times_; }
            const std::string& directory() const { return directory_; }
    };

    /// @brief Writes the k-resolved populations in the bloch gauge, @f$ P_{n}(\textbf{k}) = \rho^{bloch}_{nn}(\textbf{k}) @f$,
    /// with the band energies (one line per k point of the EDUS grid). DMk_bloch__ is local, rank 0 writes.
    void print_population_k(const std::string& filename__, const BlockMatrix<std::complex<double>>& DMk_bloch__,
                            const BandStructure& bandstructure__, const MeshGrid& kgrid__,
                            const MPIindex<3>& mpindex__, int filledbands__, const mpi::Communicator& comm__);

    /// @brief Fat bands of the populations: interpolates the density matrix on the path of kpath__,
    /// @f$ \rho(\textbf{k}) = \sum_{\textbf{R}} e^{i\textbf{k}\cdot\textbf{R}} \rho(\textbf{R}) @f$, and projects it on the bands
    /// of H__ along the same path (with the gap opened by opengap__, as in System). DMR__ is the density matrix in R space, wannier gauge (local R points).
    /// Writes filename__ (like BANDSTRUCTURE.txt, with populations instead of the wannier weights) and a gnuplot script.
    void print_fatbands_population(const std::string& filename__, const BlockMatrix<std::complex<double>>& DMR__,
                                   Operator<std::complex<double>> H__, const std::vector<std::vector<double>>& kpath__,
                                   const GridStructure& gridstructure__, int filledbands__, double opengap__,
                                   const mpi::Communicator& comm__);
} // namespace electron

#endif
