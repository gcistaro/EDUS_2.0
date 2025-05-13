
#include "mdContainers/mdContainers.hpp"
#include "MeshGrid/MeshGrid.hpp"
#include "Operator/Operator.hpp"
#include "StreamFile.hpp"
#include "ModelCoulomb/ModelCoulomb.hpp"

/// @brief Class that takes care of the calculation of the effective Hamiltonian during time propagation, with the coulomb interaction matrix elements
class Coulomb 
{
    private:
        /// Grid in R space used during the simulation
        std::shared_ptr<MeshGrid> Rgrid_;
        /// Density matrix in wannier gauge for the ground state
        Operator<std::complex<double>> DM0_;
        /// Boolean that is set to true if in the simulation we want to include the coulomb interaction in the Hamiltonian
        bool DoCoulomb_ = false;
        /// Model for the coulomb interaction, both bare and screened. We only store non-zero matrix elements
        ModelCoulomb modelcoulomb_;
        /// True if the current rank propagates R=0 term
        bool HasOrigin_;
        /// Index of the origin R=0 inside the current rank
        int index_origin_local_ = -1;
        /// The effective Hartree potential, defined as: @f[ H_{nm} = \sum_\textbf{R} V_{nm}(\textbf{R}) = 
        /// \sum_\textbf{R} \langle n\textbf{0}m\textbf{R}|V(r-r')| n\textbf{0}m\textbf{R} \rangle @f]
        mdarray<std::complex<double>, 2> Hartree;
    public:
        Coulomb(){};
        Coulomb(const int& nbnd, const std::shared_ptr<MeshGrid>& Rgrid__, const std::array<Operator<std::complex<double>>,3>& r);
        void initialize(const int& nbnd, const std::shared_ptr<MeshGrid>& Rgrid__, const std::array<Operator<std::complex<double>>,3>& r);
        void EffectiveHamiltonian(Operator<std::complex<double>>& H__, const Operator<std::complex<double>>& DM__,
                                  const bool& EraseH__);     
        void set_DM0( const Operator<std::complex<double>>& DM0__ );
        void set_DoCoulomb(const bool& DoCoulomb__);
        void set_epsilon(const double& Epsilon__);
        void set_r0(const double& r0__);
        const bool& get_DoCoulomb() const;
        bool& get_DoCoulomb();

        mdarray<std::complex<double>,3>& get_ScreenedPotential();
};



