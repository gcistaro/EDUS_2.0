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
        double r0_ = Convert(10.,Angstrom,AuLength);
        /// The macroscopic dielectric constant
        double epsilon_ = 2.;
        /// The MeshGrid of the simulation in R space
        std::shared_ptr<MeshGrid> Rgrid_;
    public:
        ModelCoulomb(){};
        ModelCoulomb(const std::array<Operator<std::complex<double>>,3>& r__, const int& dim__,
                     const std::shared_ptr<MeshGrid>& MasterRGrid__);
        void initialize(const std::array<Operator<std::complex<double>>,3>& r__, const int& dim__,
                     const std::shared_ptr<MeshGrid>& MasterRGrid_);
        void initialize_Potential( const std::shared_ptr<MeshGrid>& Rgrid__, const int& nbnd__, mdarray<std::complex<double>,3>& Potential__, 
                               const std::vector<Coordinate>& wannier_centers__, const bool& bare__ );
        std::complex<double> V(const Coordinate& r__);
        std::complex<double> W(const Coordinate& r__);
        void set_epsilon(const bool& Epsilon__);
        void set_r0(const bool& r0__);
        /// The screened interaction: ScreenedPotential_(iR, i, j)  = @f$ \langle i \textbf{0}, j \textbf{R}[iR] | W | i \textbf{0}, j \textbf{R}[iR] \rangle @f$
        mdarray<std::complex<double>,3> ScreenedPotential_;
        /// The bare interaction: BarePotential_(iR, i, j)  = @f$ \langle i \textbf{0}, j \textbf{R}[iR] | W | i \textbf{0}, j \textbf{R}[iR] \rangle @f$ 
        mdarray<std::complex<double>,3> BarePotential_;
};

#endif