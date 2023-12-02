template<class T>
Matrix<T>::Matrix(const int& nrows, const int& ncols)
{
    this->initialize(nrows, ncols);
}

template<class T>
void Matrix<T>::initialize(const int& nrows, const int& ncols)
{
    (*this).~Matrix<T>();
    assert(nrows > 0 && ncols > 0);
    this->Values.initialize({nrows, ncols});
}


//copy constructor
template<class T>
Matrix<T>::Matrix(const Matrix<T>& A){
    *this = A;
}


//copy assigment
template<class T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& m){
    (*this).~Matrix<T>();
    this->Values = m.Values;
    return *this;
}


//move constructor
template<class T>
Matrix<T>::Matrix(Matrix&& A){
    *this = A;
}


//move assignment
template<class T>
Matrix<T>& Matrix<T>::operator=(Matrix<T>&& m){
    (*this).~Matrix<T>();
    (*this).Values = m.Values;
    return *this;
}


template<class T>
const T& Matrix<T>::operator()(const int& row, const int& col) const
{
    assert ( row >= 0 && col >= 0 && row < this->get_nrows() && col < this->get_ncols() );
    return this->Values(row, col);//row major
}


template<class T>
T& Matrix<T>::operator()(const int& row, const int& col) 
{
    return const_cast<T&>(static_cast<const Matrix<T>&> (*this)(row, col));//row major
}


template<typename T>
void Matrix_gemm(Matrix<T>& OutputMatrix, const T& alpha, const Matrix<T>& InputMatrix1, const Matrix<T>& InputMatrix2, const T& beta)
{
    if(InputMatrix1.data() == nullptr || InputMatrix2.data() == nullptr){
        return;
    }
    
    //assert itis possible to do matrix multiplication
    assert( InputMatrix1.get_ncols() == InputMatrix2.get_nrows() );
    
    //initialize C if needed
    if( OutputMatrix.get_nrows() != InputMatrix1.get_nrows() ||
        OutputMatrix.get_ncols() != InputMatrix2.get_ncols())
    {
        OutputMatrix = Matrix<T>( InputMatrix1.get_nrows(), InputMatrix2.get_ncols() );
    }
    //get a notation that is consistent with cblas documentation
    auto m = InputMatrix1.get_nrows();
    auto n = InputMatrix2.get_ncols();
    auto k = InputMatrix1.get_ncols();
    //////////////////////////////////////////////////////

    gemm(m, n, k, alpha, &(InputMatrix1(0,0)), k, &(InputMatrix2(0,0)), n, beta, &(OutputMatrix(0,0)), n);
}


template<class T>
void Matrix<T>::LUdecompose(Matrix<T>& LU, lapack_int** pointer_to_ipiv) const
{
    //output: LU decomposition in a lone matrix. (upper part -> U , lower part-> L)
    //L has diagonal elements equal to 1 and are not saved; the diagonal elements are that of U.
    assert((std::is_same<T,double>::value));
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


template<typename T>
T Matrix<T>::determinant() const
{
    //in LU decomposition, we can get advantage of the fact tha the determinant of a trinagular matrix 
    //is the product of the diagonal elements, which also reprensent the eigenvalues.
    assert((*this).get_nrows() == (*this).get_ncols());
    Matrix<T> LU;
    lapack_int* ipiv;
    this->LUdecompose(LU, &ipiv);
    T determinant = 1.;
    
    for(int i=0; i < (*this).get_nrows(); ++i){
        determinant *= LU(i,i);
    }
    return determinant;
}


template<class T>
Matrix<T> Matrix<T>::inverse() const
{
    assert((*this).get_nrows() == (*this).get_ncols());
    assert((std::is_same<T,double>::value));
    assert(abs(this->determinant()) > 1.e-08);
    Matrix<T> invM;
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

template<class T>
Matrix<T> Matrix<T>::transpose() const
{
    Matrix<T> transposeM(this->get_ncols(), this->get_nrows());
    for(int irow=0; irow<this->get_ncols(); irow++){
        for(int icol=0; icol<this->get_nrows(); icol++){
            transposeM(irow, icol) = (*this)(icol, irow);
        }
    }
    return transposeM;
}

template<typename T>
void Matrix<T>::diagonalize(Matrix<T>& EigenVectors, mdarray<T,1>& EigenValues) const
{
    //note: we need to copy the matrix or we lose the info because it is overwritten with eigenvectors
    assert(this->get_nrows() == this->get_ncols());
    EigenVectors = *this;
    auto n = this->get_nrows();
    auto lda = n;
    EigenValues.initialize({this->get_nrows()});
    LAPACKE_dsyev( LAPACK_ROW_MAJOR, 'V', 'U', n, &EigenVectors(0,0), lda, &EigenValues(0) );
 
}

template<class T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& B) const
{
    Matrix<T> C;
    auto& A = static_cast<const Matrix<T>&>(*this);
    Matrix_gemm(C, 1., A, B, 0.);
    return C;
}

template<class T>
Matrix<T> Matrix<T>::operator*(T Scalar) const
{
    Matrix<T> C(this->get_nrows(), this->get_ncols());
    for(int irow=0; irow<this->get_nrows(); irow++){
        for(int icol=0; icol<this->get_ncols(); icol++){
            C(irow,icol) = Scalar*(*this)(irow,icol);
        }
    }
    return C;
}

template<class T>
Matrix<T> operator*(T Scalar, const Matrix<T> A)
{
    return A*Scalar;
}

template<class T>
Vector<T> Matrix<T>::operator*(const Vector<T>& v) const
{
    assert(v.get_NumberOfElements() == (*this).get_ncols());
    Vector<T> result(v.get_NumberOfElements());
    int m = this->get_nrows();
    int k = this->get_ncols();
    int n = 1;
    gemm(m, n, k, 1., &(*this)(0,0), k, &v(0), n, 0., &result(0), n);
    return result;

}




template<class T>
int Matrix<T>::get_nrows() const
{
    return this->Values.get_Size(0);
}


template<class T>
int Matrix<T>::get_ncols() const
{
    return this->Values.get_Size(1);
}


//destructor
template<typename T>
Matrix<T>::~Matrix()
{
    (*this).Values.~mdarray<T,2>();
}


//overloading writing matrix
/*template<class T>
std::ostream& operator<<(std::ostream& os, Matrix<T> m)
{
    for(int irow=0; irow<m.get_nrows(); irow++){
        for(int icol=0; icol<m.get_ncols(); icol++){
            os << std::setprecision(5) << std::setw(10) << m(irow,icol); 
        }
        os << std::endl;
    }
    return os;
}*/


template<class T>
std::ostream& operator<<(std::ostream& os, const Matrix<T>& m)
{
    for(int irow=0; irow<m.get_nrows(); irow++){
        for(int icol=0; icol<m.get_ncols(); icol++){
            os << std::setprecision(5) << std::setw(15) << m(irow,icol); 
        }
        os << std::endl;
    }
    return os;
}
