#include "Geometry/Matrix.hpp"

int main()
{
    Matrix<double> A(2,4);

    A(0,0) = 1.;       A(0,1) = 0.;       A(0,2) = 1.;        A(0,3) = 0.;    
    A(1,0) = 0.;       A(1,1) = 1.;       A(1,2) = 0.;        A(1,3) = 1.;

    Matrix<double> invA_(4,2);
    invA_(0,0) = 0.5;       invA_(0,1) = 0.0;          
    invA_(1,0) = 0.0;       invA_(1,1) = 0.5;       
    invA_(2,0) = 0.5;       invA_(2,1) = 0.0;          
    invA_(3,0) = 0.0;       invA_(3,1) = 0.5;       


    auto invA = A.pseudoinv();


    assert( (invA-invA_).norm() < 1.e-08 );


}