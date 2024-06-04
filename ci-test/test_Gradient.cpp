#include "Constants.hpp"
#include "Geometry/Coordinate.hpp"
#include "core/print_timing.hpp"
#include "MeshGrid/MeshGrid.hpp"
#include "kGradient/kGradient.hpp"

int main()
{   
    auto a = 2.5;      
    Matrix<double> B(3,3);

    std::cout << pi << " " << a << std::endl;
    B(0,0) = 4.*pi/(std::sqrt(3)*a)*std::cos(pi/3.);             B(0,1) = 4.*pi/(std::sqrt(3)*a)*std::cos(-pi/3.);           B(0,2) = 0.;    
    B(1,0) = 4.*pi/(std::sqrt(3)*a)*std::sin(pi/3.);             B(1,1) = 4.*pi/(std::sqrt(3)*a)*std::sin(-pi/3.);           B(1,2) = 0.;
    B(2,0) = 0.;                                                 B(2,1) = 0.;                                                B(2,2) = 10.;

    std::cout << B;
    Basis Bbasis(B);
    std::cout << "1\n";

    Coordinate::add_Basis( Bbasis, LatticeVectors(k) );
    std::cout << "2\n";
    MeshGrid mesh(k, std::array<int,3>({20, 20, 1}));
    std::cout << "3\n";

    kGradient gradient(mesh);
    std::cout << "4\n";
}
