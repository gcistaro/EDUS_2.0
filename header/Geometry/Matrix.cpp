#include <complex>
#ifndef MKL_Complex16
    #define MKL_Complex16 std::complex<double>
#endif 

#include "mkl.h"
#include "Geometry/Matrix.hpp"
template<>
void Matrix<std::complex<double>>::diagonalize(Matrix<std::complex<double>>& EigenVectors, mdarray<double,1>& EigenValues) const
{
    //note: we need to copy the matrix or we lose the info because it is overwritten with eigenvectors
    assert(this->get_nrows() == this->get_ncols());
    if( (EigenVectors.get_nrows() != this->get_nrows()) || (EigenVectors.get_ncols() != this->get_ncols()) ){
        EigenVectors.initialize(this->get_nrows(), this->get_ncols());
    }
    EigenVectors = *this;
    auto n = this->get_nrows();
    auto lda = n;
    EigenValues.initialize({this->get_nrows()});
    //LAPACKE_dsyev( LAPACK_ROW_MAJOR, 'V', 'U', n, &EigenVectors(0,0), lda, &EigenValues(0) );
    LAPACKE_zheev( LAPACK_ROW_MAJOR, 'V', 'U', n, &EigenVectors(0,0), lda, &EigenValues(0) );
}



template<>
Matrix<double> Matrix<double>::inverse() const
{
    assert((*this).get_nrows() == (*this).get_ncols());
    //assert((std::is_same<T,double>::value));
    assert(abs(this->determinant()) > 1.e-08);
    Matrix<double> invM;
    lapack_int* ipiv;
    LUdecompose(invM, &ipiv);
    //inverse
    lapack_int n = invM.get_ncols();
    lapack_int lda = n;
    LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, &invM(0,0),
                    lda, ipiv);
    delete[] ipiv;
    return invM;
}