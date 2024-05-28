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

template<>
void Matrix<double>::svd(Matrix<double>& u, Matrix<double>& vt, mdarray<double,1>& s) 
{
    //computes A = U*s*VT
    assert((*this).get_nrows() != (*this).get_ncols());

    auto jobu  = 'A'; //compute all U values
    auto jobvt = 'A'; //compute all Vt values

    lapack_int m = (*this).get_nrows();
    lapack_int n = (*this).get_ncols();
    lapack_int lda = n;
    lapack_int ldu = m;
    lapack_int ldvt = n;

    s.initialize({size_t(std::min(m,n))});          //vector with pseudo-eigenvalues
    u.initialize(ldu, m);               //left eigenvectors
    vt.initialize(ldvt, n);             //right eigenvectors (already transpose)
    mdarray<double,1> superb({size_t(std::min(m,n)-1)});
    auto info = LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'A', 'A', m, n, &((*this)(0,0)),
                               lda, s.begin().data(), &(u(0,0)), ldu,
                               &(vt(0,0)), ldvt, superb.begin().data());
    assert(info == 0);
}


template<>
Matrix<double> Matrix<double>::pseudoinv()
{
    Matrix<double> u;
    Matrix<double> vt;
    mdarray<double,1> s;
    this->svd(u, vt, s);

    //invert eigenvectors
    auto ut = u.transpose();
    auto v = vt.transpose();

    Matrix<double> pseudoinv(this->get_ncols(), this->get_nrows());
    pseudoinv.fill(0.);

    //pseudoinv = inv(vt)*1/s*inv(u)
    for( int irow = 0; irow < get_ncols(); ++irow) {
        for( int icol = 0; icol < get_nrows(); ++icol ) {
            for( int index = 0; index < s.get_TotalSize(); ++index) {
                auto sinv = ( ( std::abs(s(index)) > 1.e-08 ) ? 1./s(index) : s(index) );
                pseudoinv(irow, icol) += v(irow, index) * sinv * ut(index, icol);
            }
        }
    }
    
    return pseudoinv;
}
