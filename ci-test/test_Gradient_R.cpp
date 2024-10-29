#include "Constants.hpp"
#include "Geometry/Coordinate.hpp"
#include "core/print_timing.hpp"
#include "MeshGrid/MeshGrid.hpp"
#include "kGradient/kGradient.hpp"
#include "initialize.hpp"
#include "fftPair/fftPair.hpp"

int main()
{   
    initialize();
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
    auto A = 2.*pi*B.inverse().transpose();
    std::cout << B;
    Basis Bbasis(B);
    Basis Abasis(A);

    Coordinate::add_Basis( Bbasis, LatticeVectors(k) );
    Coordinate::add_Basis( Abasis, LatticeVectors(R) );
    MeshGrid mesh(k, std::array<int,3>({n_el, n_el, 1}));
    MeshGrid Rmesh = fftPair(mesh);
    kGradient gradient(Rmesh);

    //define function e(i 2\pi (k_0))+e(i 2\pi (k_1))
#ifdef EDUS_MPI
    mdarray<std::complex<double>, 2> Function({mesh.get_TotalSize(), 1});
    mdarray<std::complex<double>, 2> Function_R({mesh.get_TotalSize(), 1});
    mdarray<std::complex<double>, 2> DFunction_k({mesh.get_TotalSize(), 1});
    mdarray<std::complex<double>, 2> DFunction_R({mesh.get_TotalSize(), 1});
#else
    mdarray<std::complex<double>, 2> Function({1, mesh.get_TotalSize()});
    mdarray<std::complex<double>, 2> Function_R({1, mesh.get_TotalSize()});
    mdarray<std::complex<double>, 2> DFunction_k({1, mesh.get_TotalSize()});
    mdarray<std::complex<double>, 2> DFunction_R({1, mesh.get_TotalSize()});
#endif
    for( int ik=0; ik<Function.get_TotalSize(); ++ik ) {
        auto& cart = mesh[ik].get("Cartesian");
#ifdef EDUS_MPI
        Function(ik, 0) =
#else
        Function(0, ik) =
#endif
                 std::exp( 2.*pi*im*(Binv(0,0)*cart[0]+Binv(0,1)*cart[1]) )+ 
                 std::exp( 2.*pi*im*(Binv(1,0)*cart[0]+Binv(1,1)*cart[1]) );
    }
    FourierTransform ft_;
    ft_.initialize(Function, Function_R, {n_el, n_el, 1});
    ft_.fft(-1);

    FourierTransform ft_D;
    ft_D.initialize(DFunction_k, DFunction_R, {n_el, n_el, 1});

    Coordinate direction(1,0,0, "Cartesian");
    gradient.Calculate(1., DFunction_R, Function_R, direction, true);
    ft_D.fft(+1);
    for( int ik=0; ik<Function.get_TotalSize(); ++ik ) {
        auto& cart = mesh[ik].get("Cartesian");     
#ifdef EDUS_MPI
        auto& Dfunction = DFunction_k(ik, 0); 
#else
        auto& Dfunction = DFunction_k(0,ik);
#endif
        auto Analytical =  2.*pi*im*Binv(0,0)*std::exp( 2.*pi*im*(Binv(0,0)*cart[0]+Binv(0,1)*cart[1]) )
                                             +2.*pi*im*Binv(1,0)*std::exp( 2.*pi*im*(Binv(1,0)*cart[0]+Binv(1,1)*cart[1]) );
        std::cout << Dfunction << " " <<  Analytical << " " << std::abs(Dfunction-Analytical)/std::abs(Dfunction) << std::endl;
        if((std::abs(Dfunction) > 1.e-06) && std::abs(Dfunction-Analytical)/std::abs(Dfunction)*100 > 1.) {
            exit(1);
        }
    }
}
