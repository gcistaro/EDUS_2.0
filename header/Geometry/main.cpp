#include "Geometry.hpp"
#include <mkl.h>




int main()
{

    Matrix<double> M(3,3);
    M(0,0) = 1.;   M(0,1) = 0;   M(0,2) = 0;
    M(1,0) = 0;   M(1,1) = 1.;   M(1,2) = 0;
    M(2,0) = 0;   M(2,1) = 0;   M(2,2) = 1.;
    std::cout <<"M:\n"<< M << std::endl;

    Basis Cartesian;
    Cartesian.initialize(M);
    std::cout << " CAARTESIAN:: \n" << Cartesian;
    Vector<R>::add_Basis(Cartesian, "Cartesian");
    M(0,0) = +1.;   M(0,1) = +5.;   M(0,2) = +1.;
    M(1,0) = -10.;   M(1,1) = +7.;   M(1,2) = -1.;
    M(2,0) = +5.;   M(2,1) = +4.;   M(2,2) = -1.;
    Basis LatticeVectors;
    LatticeVectors.initialize(M);
    Vector<R>::add_Basis(LatticeVectors, "LatticeVectors");
    Vector<R> k1;
    k1.initialize(0.,1.,0.,"LatticeVectors");
    std::cout << "k1::: "<< k1;
}
