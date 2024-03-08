template<typename T, Space space>
void BlockMatrix<T,space>::initialize(mdarray<T,3>& Values_)//, MeshGrid& meshgrid_)
{
    Values = std::move(Values_);
    //meshgrid = std::make_shared<MeshGrid>(meshgrid_);
    initialize_submatrix();
}

template<typename T, Space space>
void BlockMatrix<T,space>::initialize_submatrix()
{
        submatrix.resize(Values.get_Size(0));
        for(auto iblock=0; iblock<Values.get_Size(0); iblock++){
            submatrix[iblock] =Matrix<T>(&(Values(iblock,0,0)),{Values.get_Size(1), Values.get_Size(2)});    
        }
}


template<typename T, Space space>
void BlockMatrix<T,space>::initialize(const size_t& nblocks, const size_t& nrows, const size_t& ncols)
{
    Values.initialize({nblocks, nrows, ncols});
    initialize_submatrix();
}

template<typename T, Space space>
void BlockMatrix<T,space>::test_submatrix()
{
    for(int ib=0; ib<this->get_nblocks(); ib++){
        for(int ir=0; ir<this->get_nrows(); ir++){
            for(int ic=0; ic<this->get_ncols(); ic++){
                std::cout << std::setw(20) << &((*this)(ib,ir,ic)); 
                std::cout << std::setw(20) << &(this->submatrix[ib](ir,ic));
                std::cout << std::setw(20) << &((*this)(ib,ir,ic)) - &(this->submatrix[ib](ir,ic));
                std::cout << std::endl;
            }
        }
    }
}


template<typename T, Space space>
BlockMatrix<T, space>::BlockMatrix(const BlockMatrix<T, space>& A)
{
    *this = A;
}

template<typename T, Space space>
const T& BlockMatrix<T, space>::operator()(const int& iblock, const int& n, const int& m) const
{
    if(iblock == -1){
        return NullValue;
    }
    return this->Values(iblock, n, m);
}

template<typename T, Space space>
T& BlockMatrix<T, space>::operator()(const int& iblock, const int& n, const int& m)
{
    return (const_cast<T&>(static_cast<BlockMatrix<T, space> const&>(*this)(iblock, n, m)));
}

template<typename T, Space space>
BlockMatrix<T, space>& BlockMatrix<T, space>::operator=(const BlockMatrix<T, space>& A)
{
    Values = A.Values;
    initialize_submatrix();
    meshgrid = A.meshgrid;
    return *this;
}

template<typename T, Space space>
BlockMatrix<T, space>::BlockMatrix(BlockMatrix<T, space>&& A)
{
    Values = std::move(A.Values);
    initialize_submatrix();
    meshgrid = std::move(A.meshgrid);
}


template<typename T, Space space>
BlockMatrix<T, space>& BlockMatrix<T, space>::operator=(BlockMatrix<T, space>&& A)
{
    Values = std::move(A.Values);
    initialize_submatrix();
    meshgrid = std::move(A.meshgrid);
    return *this;
}

template<typename T, Space space>
void BlockMatrix<T,space>::fill(const T& Scalar)
{
    PROFILE("BlockMatrix::fill");
    std::fill(this->Values.begin(), this->Values.end(), Scalar);
}

template<typename T, Space space>
Matrix<T>& BlockMatrix<T,space>::operator[](const int& iblock)
{
    return const_cast<Matrix<T>&>(static_cast<const BlockMatrix<T,space>&>(*this)[iblock]);
}

template<typename T, Space space>
const Matrix<T>& BlockMatrix<T,space>::operator[](const int& iblock) const
{
    if(iblock < 0){
        return EmptyMatrix;
    }
    return submatrix[iblock];
}

template<typename T, Space space>
Matrix<T>& BlockMatrix<T,space>::operator[](const Coordinate<space>& Point)
{
    return const_cast<Matrix<T>&>(static_cast<const BlockMatrix<T,space>&>(*this)[Point]);
}

template<typename T, Space space>
const Matrix<T>& BlockMatrix<T,space>::operator[](const Coordinate<space>& Point) const
{
    int iblock = (*meshgrid).find(Point);
    return (*this)[iblock];
}


template<typename T, Space space>
void multiply(BlockMatrix<T,space>& Output, T Scalar, const BlockMatrix<T,space>& Input1, const BlockMatrix<T,space>& Input2 )
{
    for(int iblock=0; iblock<Output.get_nblocks(); iblock++){
        Matrix_gemm(Output[iblock], Scalar, Input1[iblock], Input2[iblock], T(0.));
    }
}

template<typename T, Space space>
void multiply(BlockMatrix<T,space>& Output, T Scalar, const BlockMatrix<T,space>& Input1, const BlockMatrix<T,space>& Input2, T Scalar2 )
{
    for(int iblock=0; iblock<Output.get_nblocks(); iblock++){
        Matrix_gemm(Output[iblock], Scalar, Input1[iblock], Input2[iblock], Scalar2);
    }
}


template<typename T, Space space, typename U>
void convolution1(BlockMatrix<T,space>& Output, U Scalar, const BlockMatrix<T,space>& Input1, const BlockMatrix<T,space>& Input2 )
{
    assert(Output.get_nblocks()!=0);

    //Output.fill(0.);
    auto& ci = MeshGrid<space>::ConvolutionIndex1[{Output.meshgrid->get_id(), Input1.meshgrid->get_id(), Input2.meshgrid->get_id()}];

    if(ci.get_Size(0) == 0 ){
        MeshGrid<space>::Calculate_ConvolutionIndex1(*(Output.meshgrid), *(Input1.meshgrid), *(Input2.meshgrid));

    }
    PROFILE_START("BlockMatrix::convolution");
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

template<typename T, Space space, typename U>
void convolution2(BlockMatrix<T,space>& Output, U Scalar, const BlockMatrix<T,space>& Input1, const BlockMatrix<T,space>& Input2 )
{
    assert(Output.get_nblocks()!=0);

    //Output.fill(0.);
    auto& ci = MeshGrid<space>::ConvolutionIndex2[{Output.meshgrid->get_id(), Input1.meshgrid->get_id(), Input2.meshgrid->get_id()}];

    if(ci.get_Size(0) == 0 ){
        MeshGrid<space>::Calculate_ConvolutionIndex2(*(Output.meshgrid), *(Input1.meshgrid), *(Input2.meshgrid));

    }
    PROFILE_START("BlockMatrix::convolution");
    for(int iblock_o=0; iblock_o<Output.get_nblocks(); iblock_o++){
        for(int iblock_i1=0; iblock_i1<Input1.get_nblocks(); iblock_i1++){
            auto& iblock_i2 = ci(iblock_o, iblock_i1);
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



template<typename T, Space space>
BlockMatrix<T,space>::~BlockMatrix()
{
    Values.~mdarray<T,3>();
}


template<typename T, Space space>
void BlockMatrix<T,space>::diagonalize(std::vector<mdarray<double,1>>& Eigenvalues,
                 BlockMatrix<std::complex<double>,space>& Eigenvectors) const
{
    //assert(space_ == k);
    Eigenvectors.initialize(this->get_nblocks(), this->get_nrows(), this->get_ncols());
    std::cout << "eigenvectors test : " << &(Eigenvectors[0](0,0)) << " "<< &(Eigenvectors(0,0,0)) << std::endl;
    Eigenvalues.resize(this->get_nblocks());
    for(int ik=0; ik<this->get_nblocks(); ik++){
        Eigenvalues[ik].initialize({this->get_nrows()});
        ((*this)[ik]).diagonalize(Eigenvectors[ik], Eigenvalues[ik]);
    }
}


template<typename T, Space space>
auto max(const BlockMatrix<T,space>& m)
{
    //double max = 0;
//    for(int iblock=0; iblock<m.get_nblocks(); ++iblock){
//        for(int irow=0; irow<m.get_nrows(); ++irow){
//            for(int icol=0; icol<m.get_ncols(); ++icol){
//                //HR(ci(iblock, 0), irow, icol) = H0R(iblock, irow, icol);
//                if(std::abs(m(iblock, irow, icol)) > max){
//                    max =std::abs( m(iblock, irow, icol));
//                }
//            }
//        }
//    }            
    return std::max_element(m.begin(), m.end(), [&](const T& a1, const T& a2){ return std::abs(a1) < std::abs(a2);});
}

template<class T, Space space>
std::ostream& operator<<(std::ostream& os, const BlockMatrix<T, space>& m)
{
    for(int i=0; i<m.get_nblocks(); i++){
        std::cout << (*((const_cast<BlockMatrix<T,space>&>(m)).get_MeshGrid()))[i].get("LatticeVectors") << m[i] << std::endl;
    }
    return os;
}
