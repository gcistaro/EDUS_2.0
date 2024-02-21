#include "Geometry/Matrix.hpp"

template<typename T, typename T_>
Matrix<T> diagonal(mdarray<T_, 1> DiagonalElements)
{
    Matrix<T> Output(DiagonalElements.get_Size(0), DiagonalElements.get_Size(0));
    Output.fill(0);
    for(int el=0; el<DiagonalElements.get_Size(0); el++){
        Output(el, el) = DiagonalElements(el);
    }
    return Output;
}


int main()
{

    Matrix<std::complex<double>> M(3,3);

    M(0,0) = +1.;   M(0,1) = -10.;   M(0,2) = +5.;
    M(1,0) = -10.;   M(1,1) = +7.;   M(1,2) = +4.;
    M(2,0) = +5.;   M(2,1) = +4.;   M(2,2) = -1.;


    Matrix<std::complex<double>> U; 
    mdarray<double,1> E;
    M.diagonalize(U,E); 
    Matrix<std::complex<double>> Lambda = diagonal<std::complex<double>>(E);
    auto Lambda_ = U.transpose()*M*U;
    std::cout << "Ut.M.U  \n" << Lambda_ << std::endl;
    std::cout << "Lambda:  \n" << Lambda << std::endl;

    if( (Lambda-Lambda_).norm() > 1.e-08){
        exit(1);
    } 

}