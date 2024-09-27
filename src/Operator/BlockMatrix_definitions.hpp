template<typename T>
Matrix<T> BlockMatrix<T>::EmptyMatrix = Matrix<T>();

template<typename T>
void BlockMatrix<T>::initialize(const Space& space, mdarray<T,3>& Values_)//, MeshGrid& meshgrid_)
{
    Values = std::move(Values_);
    //meshgrid = std::make_shared<MeshGrid>(meshgrid_);
    initialize_submatrix();
}

template<typename T>
void BlockMatrix<T>::initialize_submatrix()
{
        submatrix.resize(Values.get_Size(0));

        for(auto iblock=0; iblock<Values.get_Size(0); iblock++){
            submatrix[iblock] =Matrix<T>(&(Values(iblock,0,0)),{Values.get_Size(1), Values.get_Size(2)});    
        }
}


template<typename T>
void BlockMatrix<T>::initialize(const Space& space__, const int& nblocks, const int& nrows, const int& ncols)
{
    Values.initialize({nblocks, nrows, ncols});
    space = space__;
    initialize_submatrix();
}

template<typename T>
void BlockMatrix<T>::test_submatrix()
{
    for(int ib=0; ib<this->get_nblocks(); ib++){
        assert( (*this)[ib].get_nrows() == this->get_nrows() );
        assert( (*this)[ib].get_ncols() == this->get_ncols() );

        for(int ir=0; ir<this->get_nrows(); ir++){
            for(int ic=0; ic<this->get_ncols(); ic++){
                assert( &((*this)(ib,ir,ic)) == &(this->submatrix[ib](ir,ic)) );
            }
        }
    }
}


template<typename T>
BlockMatrix<T>::BlockMatrix(const BlockMatrix<T>& A)
{
    *this = A;
}

template<typename T>
inline const T& BlockMatrix<T>::operator()(const int& iblock, const int& n, const int& m) const
{
    if(iblock == -1){
        return NullValue;
    }
    return this->Values(iblock, n, m);
}

template<typename T>
inline T& BlockMatrix<T>::operator()(const int& iblock, const int& n, const int& m)
{
    return (const_cast<T&>(static_cast<BlockMatrix<T> const&>(*this)(iblock, n, m)));
}

template<typename T>
inline T& BlockMatrix<T>::operator()(const int& i)
{
    return this->Values[i];
}

template<typename T>
inline const T& BlockMatrix<T>::operator()(const int& i) const
{
    return (const_cast<T&>(static_cast<BlockMatrix<T> const&>(*this)(i)));
}

template<typename T>
BlockMatrix<T>& BlockMatrix<T>::operator=(const BlockMatrix<T>& A)
{
    Values = A.Values;
    space = A.space;
    initialize_submatrix();
    meshgrid = A.meshgrid;
    return *this;
}

template<typename T>
BlockMatrix<T>::BlockMatrix(BlockMatrix<T>&& A)
{
    Values = std::move(A.Values);
    space = std::move(A.space);
    initialize_submatrix();
    meshgrid = std::move(A.meshgrid);
}


template<typename T>
BlockMatrix<T>& BlockMatrix<T>::operator=(BlockMatrix<T>&& A)
{
    Values = std::move(A.Values);
    space = std::move(A.space);
    initialize_submatrix();
    meshgrid = std::move(A.meshgrid);
    return *this;
}

template<typename T>
void BlockMatrix<T>::fill(const T& Scalar)
{
    PROFILE("BlockMatrix::fill");
    std::fill(this->Values.begin(), this->Values.end(), Scalar);
}

template<typename T>
Matrix<T>& BlockMatrix<T>::operator[](const int& iblock)
{
    return const_cast<Matrix<T>&>(static_cast<const BlockMatrix<T>&>(*this)[iblock]);
}

template<typename T>
const Matrix<T>& BlockMatrix<T>::operator[](const int& iblock) const
{
    if(iblock < 0){
        return EmptyMatrix;
    }
    return submatrix[iblock];
}

template<typename T>
Matrix<T>& BlockMatrix<T>::operator[](const Coordinate& Point)
{
    return const_cast<Matrix<T>&>(static_cast<const BlockMatrix<T>&>(*this)[Point]);
}

template<typename T>
const Matrix<T>& BlockMatrix<T>::operator[](const Coordinate& Point) const
{
    int iblock = (*meshgrid).find(Point);
    return (*this)[iblock];
}


template<typename T>
void multiply(BlockMatrix<T>& Output, T Scalar, const BlockMatrix<T>& Input1, const BlockMatrix<T>& Input2 )
{
    multiply(Output, Scalar, Input1, Input2, T(0.));
}

template<typename T>
void multiply(BlockMatrix<T>& Output, T Scalar, const BlockMatrix<T>& Input1, const BlockMatrix<T>& Input2, T Scalar2 )
{
    assert(Output.get_nblocks() == Input1.get_nblocks() && Input1.get_nblocks() == Input2.get_nblocks());
    #pragma omp parallel for schedule(static)
    for(int iblock=0; iblock<Output.get_nblocks(); iblock++){
        Matrix_gemm(Output[iblock], Scalar, Input1[iblock], Input2[iblock], Scalar2);
    }
}


template<typename T, typename U>
void convolution(BlockMatrix<T>& Output, U Scalar, const BlockMatrix<T>& Input1, const BlockMatrix<T>& Input2 )
{
    assert(Output.get_nblocks()!=0);

    //Output.fill(0.);
    auto& ci = MeshGrid::ConvolutionIndex[{Output.meshgrid->get_id(), Input1.meshgrid->get_id(), Input2.meshgrid->get_id()}];

    if(ci.get_Size(0) == 0 ){
        MeshGrid::Calculate_ConvolutionIndex(*(Output.meshgrid), *(Input1.meshgrid), *(Input2.meshgrid));

    }
    PROFILE_START("BlockMatrix::convolution");

    #pragma omp parallel for schedule(static)
    for(int iblock_o=0; iblock_o<Output.get_nblocks(); iblock_o++){
        for(int iblock_i2=0; iblock_i2<Input2.get_nblocks(); iblock_i2++){
            auto& iblock_i1 = ci(iblock_o, iblock_i2);
            Matrix_gemm(Output[iblock_o], Scalar+im*0., Input1[iblock_i1], Input2[iblock_i2], 1.+im*0.);
        }
    }

    //auto maximum1 = max(Input1);
    //auto maximum2 = max(Input2);
    //auto maximum3 = max(Output);
    //std::cout << "MAX INPUT1:  " << maximum1 << std::endl;
    //std::cout << "MAX INPUT2:  " << maximum2 << std::endl;
    //std::cout << "MAX OUTPUT:  " << maximum3 << std::endl;
    PROFILE_STOP("BlockMatrix::convolution");
}

template<typename T_, typename U>
void commutator(BlockMatrix<T_>& Output, U Scalar, const BlockMatrix<T_>& Input1, const BlockMatrix<T_>& Input2, const bool& Erase_Output = true)
{
    PROFILE("Commutator");
    assert( Output.get_space() == Input1.get_space() );
    assert( Input2.get_space() == Input1.get_space() );

    switch( Output.get_space() ){
        case(R):
        {
            convolution(Output, Scalar, Input1, Input2);
            convolution(Output, -Scalar, Input2, Input1);
            break;
        }
        case(k):
        {
            multiply(Output, Scalar, Input1, Input2, double(!Erase_Output) + im*0.);
            multiply(Output, -Scalar, Input2, Input1, 1.+im*0.);
            break;
        }
    }
}



template<typename T>
BlockMatrix<T>::~BlockMatrix()
{
    Values.~mdarray<T,3>();
}


template<typename T>
void BlockMatrix<T>::diagonalize(std::vector<mdarray<double,1>>& Eigenvalues,
                 BlockMatrix<std::complex<double>>& Eigenvectors) const
{
    //assert(space_ == k);
    Eigenvectors.initialize(k, this->get_nblocks(), this->get_nrows(), this->get_ncols());
    Eigenvalues.resize(this->get_nblocks());
    for(int ik=0; ik<this->get_nblocks(); ik++){
        Eigenvalues[ik].initialize({this->get_nrows()});
        ((*this)[ik]).diagonalize(Eigenvectors[ik], Eigenvalues[ik]);
    }
}

template<typename T>
bool BlockMatrix<T>::is_hermitian()
{
    assert(space == k);
    bool hermitian = true;
    for(int ik=0; ik<this->get_nblocks(); ++ik) {
        for( int irow=0; irow<this->get_nrows(); irow++ ) {
            for(int icol=irow; icol < this->get_nrows(); ++icol ) {
                //std::cout << (*this)[ik](irow, icol) << "    " << (*this)[ik](icol, irow) << "  " ;
                //std::cout << std::abs( (*this)[ik](irow, icol) - std::conj((*this)[ik](icol, irow)) ) << std::endl;  
                hermitian = hermitian && ( std::abs( (*this)[ik](irow, icol) - std::conj((*this)[ik](icol, irow)) ) < 1.e-15) ;
            }
        }
    }
    //std::cout <<  ( hermitian ? "IS HERMITIAN!!":"IS NOT HERMITIAN :(" ) << std::endl;
    return hermitian;
}

template<typename T>
void BlockMatrix<T>::make_hermitian()
{
    switch(space) 
    {
        case(Space::k) :
        {
            #pragma omp parallel for
            for(int ik=0; ik<this->get_nblocks(); ++ik) {
                for( int irow=0; irow<this->get_nrows(); irow++ ) {
                    for(int icol=irow; icol < this->get_nrows(); ++icol ) {
                        //std::cout << (*this)[ik](irow, icol) << "    " << (*this)[ik](icol, irow) << "  " ;
                        //std::cout << std::abs( (*this)[ik](irow, icol) - std::conj((*this)[ik](icol, irow)) ) << std::endl;
                        auto term = ((*this)(ik,irow,icol) + std::conj((*this)(ik,icol,irow)))/2.;
                        (*this)(ik,irow,icol) = term;
                        (*this)(ik,icol,irow) = std::conj(term);
                    }
                }
            }
            break;
        }
        case(Space::R) :
        {
            std::cout << "Implement me\n";
            exit(1);
        }
        
    }

}


template<typename T>
void BlockMatrix<T>::make_antihermitian()
{
    #pragma omp parallel for
    for(int ik=0; ik<this->get_nblocks(); ++ik) {
        for( int irow=0; irow<this->get_nrows(); irow++ ) {
            for(int icol=irow; icol < this->get_nrows(); ++icol ) {
                //std::cout << (*this)[ik](irow, icol) << "    " << (*this)[ik](icol, irow) << "  " ;
                //std::cout << std::abs( (*this)[ik](irow, icol) - std::conj((*this)[ik](icol, irow)) ) << std::endl;
                auto term = ((*this)(ik,irow,icol) - std::conj((*this)(ik,icol,irow)))/2.;
                (*this)(ik,irow,icol) = term;
                (*this)(ik,icol,irow) = -std::conj(term);
            }
        }
    }

}


template<typename T>
void BlockMatrix<T>::cut(const double& threshold__)
{
    #pragma omp parallel for
    for(int ik=0; ik<this->get_nblocks(); ++ik) {
        for( int irow=0; irow<this->get_nrows(); irow++ ) {
            for(int icol=irow; icol < this->get_nrows(); ++icol ) {
                (*this)(ik,irow,icol) = ( std::abs( (*this)(ik,irow,icol) ) > threshold ? (*this)(ik,irow,icol) : 0.);
            }
        }
    }
}

template<typename T>
auto max(const BlockMatrix<T>& m)
{
    return std::max_element(m.begin(), m.end(), [&](const T& a1, const T& a2){ return std::abs(a1) < std::abs(a2);});
}

template<class T>
std::ostream& operator<<(std::ostream& os, const BlockMatrix<T>& m)
{
    std::cout << m.get_nblocks();
    for(int i=0; i<m.get_nblocks(); i++){
        os << i << " " << (*((const_cast<BlockMatrix<T>&>(m)).get_MeshGrid()))[i].get(LatticeVectors(m.get_space())) << m[i] << std::endl;
    }
    return os;
}


