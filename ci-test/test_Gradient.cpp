#include "Constants.hpp"
#include "Geometry/Coordinate.hpp"
#include "core/print_timing.hpp"
#include "MeshGrid/MeshGrid.hpp"
#include "kGradient/kGradient.hpp"

int main()
{   
    int n_el = 100;
    auto a = 2.5;      
    Matrix<double> B(3,3);
/*
    B(0,0) = 1.;             B(0,1) = 0.;           B(0,2) = 0.;    
    B(1,0) = 0.;             B(1,1) = 1.;           B(1,2) = 0.;
    B(2,0) = 0.;             B(2,1) = 0.;           B(2,2) = 10.;
*/
    B(0,0) = 4.*pi/(std::sqrt(3)*a)*std::cos(pi/3.);             B(0,1) = 4.*pi/(std::sqrt(3)*a)*std::cos(-pi/3.);           B(0,2) = 0.;    
    B(1,0) = 4.*pi/(std::sqrt(3)*a)*std::sin(pi/3.);             B(1,1) = 4.*pi/(std::sqrt(3)*a)*std::sin(-pi/3.);           B(1,2) = 0.;
    B(2,0) = 0.;                                                 B(2,1) = 0.;                                                B(2,2) = 10.;

    auto Binv = B.inverse();
    std::cout << B;
    Basis Bbasis(B);

    Coordinate::add_Basis( Bbasis, LatticeVectors(k) );
    MeshGrid mesh(k, std::array<int,3>({n_el, n_el, 1}));
    kGradient gradient(mesh);

    //define function e(i 2\pi (k_0))+e(i 2\pi (k_1))
    mdarray<std::complex<double>, 1> Function({mesh.get_TotalSize()});
    for( int ik=0; ik<Function.get_TotalSize(); ++ik ) {
        auto& cart = mesh[ik].get("Cartesian");
        Function(ik) = std::exp( 2.*pi*im*(Binv(0,0)*cart[0]+Binv(0,1)*cart[1]) )+ 
                       std::exp( 2.*pi*im*(Binv(1,0)*cart[0]+Binv(1,1)*cart[1]) );
    }
    mdarray<std::complex<double>, 1> Derivative({mesh.get_TotalSize()});
    Coordinate direction(1,0,0, "Cartesian");
    gradient.Calculate(Derivative, Function, direction, true);
    for( int ik=0; ik<Function.get_TotalSize(); ++ik ) {
        auto& cart = mesh[ik].get("Cartesian");       
        auto Analytical =  2.*pi*im*Binv(0,0)*std::exp( 2.*pi*im*(Binv(0,0)*cart[0]+Binv(0,1)*cart[1]) )
                                             +2.*pi*im*Binv(1,0)*std::exp( 2.*pi*im*(Binv(1,0)*cart[0]+Binv(1,1)*cart[1]) );
        std::cout << Derivative(ik) << " " <<  Analytical << " " << std::abs(Derivative(ik)-Analytical)/std::abs(Derivative(ik)) << std::endl;
        if((std::abs(Derivative(ik)) > 1.e-06) && std::abs(Derivative(ik)-Analytical)/std::abs(Derivative(ik))*100 > 1.) {
            exit(1);
        }
    }
}
