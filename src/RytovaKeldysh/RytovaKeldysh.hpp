#include <math.h>
#include "Constants.hpp"
#include "ConvertUnits.hpp"
#include "Geometry/Coordinate.hpp"
#include "Operator/Operator.hpp"
#include "fftPair/fftPair.hpp"

//computes H_v(x)
double struve(const double& x, const double& v);
//computes Y_v(x)
//this is contained in math.h

class RytovaKeldysh
{
    private:
        enum Dimensionality {twoD, threeD} dim;
        double r0=Convert(10.,Angstrom,AuLength);
        double epsilon = 2.;
        double Area=1.;

        std::shared_ptr<MeshGrid> Rgrid;
    public:
        RytovaKeldysh(){};
        RytovaKeldysh(const std::array<Operator<std::complex<double>>,3>& r, const int& dim_,
                     const std::shared_ptr<MeshGrid>& MasterRGrid);
        void initialize(const std::array<Operator<std::complex<double>>,3>& r, const int& dim_,
                     const std::shared_ptr<MeshGrid>& MasterRGrid);
        //std::complex<double> W(const Coordinate& q);
        std::complex<double> Potential(const Coordinate& r);

        mdarray<std::complex<double>,3> TB;
};
