template<class T>
Matrix<T>::Matrix(const size_t& nrows, const size_t& ncols)
{
    this->initialize(nrows, ncols);
}

template<class T>
void Matrix<T>::initialize(const size_t& nrows, const size_t& ncols)
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
    if( this->get_nrows() != m.get_nrows() || this->get_ncols() != m.get_ncols() ){
        (*this).~Matrix<T>();
    }
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
Matrix<T>& Matrix<T>::operator=(Matrix<T>&& m)
{
    (*this).~Matrix<T>();
    (*this).Values = std::move(m.Values);
    return *this;
}


template<class T>
Matrix<T>::Matrix(T* Ptr, const std::array<size_t,2>& Size_)
{
    Values.initialize(Ptr, Size_);
}

template<class T>
void Matrix<T>::fill(const T& filling_constant)
{
    Values.fill(filling_constant);
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


template<typename T, typename T_>
void Matrix_gemm(Matrix<T>& OutputMatrix, const T_& alpha, const Matrix<T>& InputMatrix1, const Matrix<T>& InputMatrix2, const T_& beta)
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
    size_t m = (*this).get_nrows();
    size_t n = (*this).get_ncols();
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


template<class T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& B) const
{
    Matrix<T> C;
    auto& A = static_cast<const Matrix<T>&>(*this);
    Matrix_gemm(C, 1., A, B, 0.);
    return C;
}

template<class T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& B) const
{
    assert( this->get_nrows() == B.get_nrows() && this->get_ncols() == B.get_ncols() );
    Matrix<T> C(this->get_nrows(), this->get_ncols());
    for(int irow=0; irow<C.get_nrows(); ++irow){
        for(int icol=0; icol<C.get_ncols(); ++icol){
            C(irow, icol) = (*this)(irow, icol) + B(irow, icol); 
        }
    }
    return C;
}

template<class T>
Matrix<T> Matrix<T>::operator-() const
{
    Matrix<T> C(this->get_nrows(), this->get_ncols());
    for(int irow=0; irow<C.get_nrows(); ++irow){
        for(int icol=0; icol<C.get_ncols(); ++icol){
            C(irow, icol) = -(*this)(irow, icol);
        }
    }
    return C;
}

template<class T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& B) const
{    
    return (*this) + (-B);
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
size_t Matrix<T>::get_nrows() const
{
    return this->Values.get_Size(0);
}


template<class T>
size_t Matrix<T>::get_ncols() const
{
    return this->Values.get_Size(1);
}


//destructor
template<typename T>
Matrix<T>::~Matrix()
{
    (*this).Values.~mdarray<T,2>();
}

template<typename T>
double Matrix<T>::norm() const
{
    double Norm=0;
    for(auto& val : (*this)){
        Norm += abs(val)*abs(val);
    }
    return Norm;
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


