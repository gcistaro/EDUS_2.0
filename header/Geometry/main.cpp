#include "Vector.hpp"
#include "Matrix.hpp"
#include "Coordinate.hpp"
#include <mkl.h>




int main()
{

    Matrix<double> M(3,3);
/*

    M(0,0) = 1.;   M(0,1) = 1;   M(0,2) = 1;
    M(1,0) = 1;   M(1,1) = -1.;   M(1,2) = 1;
    M(2,0) = 1;   M(2,1) = 1;   M(2,2) = -1.;
    std::cout <<"M:\n"<< M << std::endl;
    std::cout <<"LU: \n" << M.LUdecompose() << std::endl;
    std::cout << "determinant:\n" << M.determinant() << std::endl;
    std::cout << M.inverse() << std::endl;
    std::cout << std::endl << M.inverse()*M << std::endl;
    Vector<double> v(3);
    v(0) = -2; v(1) = 0; v(2) = 1; 
    auto v1 = M*v;
    Vector<double> v2;
*/
    //v2.initialize_n(1,2,3);
    //std::cout << v;
    //std::cout << v1;
    //std::cout << v2;

    M(0,0) = +1.;   M(0,1) = -10.;   M(0,2) = +5.;
    M(1,0) = -10.;   M(1,1) = +7.;   M(1,2) = +4.;
    M(2,0) = +5.;   M(2,1) = +4.;   M(2,2) = -1.;


    Matrix<double> U; 
    mdarray<double,1> E;
    M.diagonalize(U,E); 
    std::cout << "Ut.M.U  \n" << U.transpose()*M*U << std::endl;
    std::cout << "ENergy \n " << E;    
/*
    auto invM = M.inverse();   
    std::cout << "I have done it! This is invM" << invM << std::endl;
    std::cout << "verification \n " <<invM*M << "\n" << M*invM ;

    Basis Cartesian;
    std::cout << "\n\n\nINitializing Cartesian basis\n\n\n";
    Cartesian.initialize(M);

    std::cout << " CAARTESIAN:: \n" << Cartesian;
    std::cout << Cartesian.get_M()*Cartesian.get_invM() << std::endl;
    Coordinate<R>::add_Basis(Cartesian, "Cartesian");
    M(0,0) = +1.;   M(0,1) = +5.;   M(0,2) = +1.;
    M(1,0) = -10.;   M(1,1) = +7.;   M(1,2) = -1.;
    M(2,0) = +5.;   M(2,1) = +4.;   M(2,2) = -1.;
    Basis LatticeVectors;
    LatticeVectors.initialize(M);
    Coordinate<R>::add_Basis(LatticeVectors, "LatticeVectors");
    Coordinate<R> k1;
    k1.initialize(0.,1.,0.,"LatticeVectors");
    std::cout << "k1::: "<< k1;
*/
}
