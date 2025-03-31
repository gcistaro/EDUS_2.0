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
    LAPACKE_zheevd( LAPACK_ROW_MAJOR, 'V', 'U', n, &EigenVectors(0,0), lda, &EigenValues(0) );
}

template<>
void Matrix<double>::LUdecompose(Matrix<double>& LU, lapack_int** pointer_to_ipiv) const
{
    //output: LU decomposition in a lone matrix. (upper part -> U , lower part-> L)
    //L has diagonal elements equal to 1 and are not saved; the diagonal elements are that of U.
    LU = *this;
    int m = (*this).get_nrows();
    int n = (*this).get_ncols();
    lapack_int lda = n;
    //if(*pointer_to_ipiv != nullptr){
    //	    delete[] *pointer_to_ipiv;
    //}
    *pointer_to_ipiv= new lapack_int[n];
    
    //LU decomposition
    LAPACKE_dgetrf(LAPACK_ROW_MAJOR, m, n, 
                   &(LU(0,0)), lda, *pointer_to_ipiv);  
}

template<>
void Matrix<std::complex<double>>::LUdecompose(Matrix<std::complex<double>>& LU, lapack_int** pointer_to_ipiv) const
{
    //output: LU decomposition in a lone matrix. (upper part -> U , lower part-> L)
    //L has diagonal elements equal to 1 and are not saved; the diagonal elements are that of U.
    LU = *this;
    int m = (*this).get_nrows();
    int n = (*this).get_ncols();
    lapack_int lda = n;
    //if(*pointer_to_ipiv != nullptr){
    //	    delete[] *pointer_to_ipiv;
    //}
    *pointer_to_ipiv= new lapack_int[n];
    
    //LU decomposition
    auto info = LAPACKE_zgetrf(LAPACK_ROW_MAJOR, m, n, 
                   &(LU(0,0)), lda, *pointer_to_ipiv);  
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
Matrix<std::complex<double>> Matrix<std::complex<double>>::inverse() const
{
    assert((*this).get_nrows() == (*this).get_ncols());
    //assert((std::is_same<T,double>::value));
    assert(abs(this->determinant()) > 1.e-08);
    Matrix<std::complex<double>> invM;
    lapack_int* ipiv;
    LUdecompose(invM, &ipiv);
    //inverse
    lapack_int n = invM.get_ncols();
    lapack_int lda = n;
    LAPACKE_zgetri(LAPACK_ROW_MAJOR, n, &invM(0,0),
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

    s.initialize({std::min(m,n)});          //vector with pseudo-eigenvalues
    u.initialize(ldu, m);               //left eigenvectors
    vt.initialize(ldvt, n);             //right eigenvectors (already transpose)
    mdarray<double,1> superb({std::min(m,n)-1});
    
    //copy *this to avoid overwriting
    auto A_svd = *this;        
    auto info = LAPACKE_dgesvd(LAPACK_ROW_MAJOR, jobu, jobvt, m, n, &(A_svd(0,0)),
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
    //1. A*A^{-1}*A = A
    //assert ( ( (*this)*pseudoinv*(*this) - (*this) ).norm() < 1.e-07 );
    //2. A^{-1}*A*A^{-1} = A^{-1}
    //assert ( ( pseudoinv*(*this)*pseudoinv - pseudoinv ).norm() < 1.e-07 );
    //3. (A*A^{-1})^H = A*A^{-1}
    //4. (A^{-1}*A)^H = A^{-1}*A
    
    return pseudoinv;
}
