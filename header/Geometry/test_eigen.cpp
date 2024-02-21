#include "Matrix.hpp"

int main()
{

    Matrix<double> M(3,3);

    M(0,0) = +1.;   M(0,1) = -10.;   M(0,2) = +5.;
    M(1,0) = -10.;   M(1,1) = +7.;   M(1,2) = +4.;
    M(2,0) = +5.;   M(2,1) = +4.;   M(2,2) = -1.;


    Matrix<double> U; 
    mdarray<double,1> E;
    M.diagonalize(U,E); 
    std::cout << "Ut.M.U  \n" << U.transpose()*M*U << std::endl;
    std::cout << "ENergy \n " << E;    
}
