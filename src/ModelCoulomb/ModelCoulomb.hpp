#include <math.h>
#include "Constants.hpp"
#include "ConvertUnits.hpp"
#include "Geometry/Coordinate.hpp"
#include "Operator/Operator.hpp"
#include "fftPair/fftPair.hpp"

#ifndef __MODELCOULOMB_HPP
#define __MODELCOULOMB_HPP
//computes H_v(x)
double struve(const double& x, const double& v);
//computes Y_v(x)
//this is contained in math.h

/// @brief Class for modelling matrix elements of Coulomb interaction in Point-Like approximation. We only store the non-zero elements
class ModelCoulomb
{
    private:
        /// dimensionality of the system -> twoD (2D) monolayer or threeD (3D) bulk material
        enum Dimensionality {twoD, threeD} dim_;
        /// The r0 of the RytovaKeldysh model
        std::array<double, 3> r0_;
        /// The average of the r0 parameters over the 3 dimensions
        double r0_avg_;
        /// The macroscopic dielectric constant
        double epsilon_ = 2.;
        /// The MeshGrid of the simulation in R space
        std::shared_ptr<MeshGrid> Rgrid_;
        /// minimum distance between two wannier centers, for saturations of bare and screened coulomb
        Coordinate min_distance_;
        /// Norm of the minimum distance
        double min_distance_norm_;
    public:
        ModelCoulomb(){};
        ModelCoulomb(const std::array<Operator<std::complex<double>>,3>& r__, const int& dim__,
                     const std::shared_ptr<MeshGrid>& MasterRGrid__, const bool& read_interaction__);
        void initialize(const std::array<Operator<std::complex<double>>,3>& r__, const int& dim__,
                     const std::shared_ptr<MeshGrid>& MasterRGrid_, const bool& read_interaction__);
        void initialize_Potential( const std::shared_ptr<MeshGrid>& Rgrid__, const int& nbnd__, mdarray<std::complex<double>,3>& Potential__, 
                               const std::vector<Coordinate>& wannier_centers__, const bool& bare__ );
        void initialize_Potential(const std::shared_ptr<MeshGrid>& Rgrid__, const int& nbnd__, mdarray<std::complex<double>,3>& Potential__,
                                const bool& bare);
        std::complex<double> V(const Coordinate& r__);
        std::complex<double> W(const Coordinate& r__);
        void set_epsilon(const double& Epsilon__);
        void set_r0(const std::vector<double>& r0__);
        std::array<double, 3>& get_r0();
        double get_r0_avg();
        /// The screened interaction: ScreenedPotential_(iR, i, j)  = @f$ \langle i \textbf{0}, j \textbf{R}[iR] | W | i \textbf{0}, j \textbf{R}[iR] \rangle @f$
        mdarray<std::complex<double>,3> ScreenedPotential_;
        /// The bare interaction: BarePotential_(iR, i, j)  = @f$ \langle i \textbf{0}, j \textbf{R}[iR] | W | i \textbf{0}, j \textbf{R}[iR] \rangle @f$ 
        mdarray<std::complex<double>,3> BarePotential_;

        mdarray<std::complex<double>,3>& get_ScreenedPotential();
};

#endif