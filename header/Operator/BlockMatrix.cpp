template<typename T>
void BlockMatrix<T>::initialize(mdarray<T,3>& Values_)//, MeshGrid& meshgrid_)
{
    Values = std::move(Values_);
    //meshgrid = std::make_shared<MeshGrid>(meshgrid_);

    submatrix.resize(Values.get_Size(0));
    auto Size1 = std::array<size_t,2>{Values.get_Size(1), Values.get_Size(2)};
    for(auto iblock=0; iblock<Values.get_Size(0); iblock++){
        submatrix[iblock] =Matrix<T>(mdarray<T,2>(&(Values(iblock,0,0)),{Values.get_Size(1), Values.get_Size(2)}));
    }
}

template<typename T>
void BlockMatrix<T>::fill(const T& Scalar)
{
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
void multiply(BlockMatrix<T>& Output, T Scalar, const BlockMatrix<T>& Input1, const BlockMatrix<T>& Input2 )
{
    for(int iblock=0; iblock<Output.get_nblocks(); iblock++){
        Matrix_gemm(Output[iblock], Scalar, Input1[iblock], Input2[iblock], 0.);
    }
}

template<typename T, typename U>
void convolution(BlockMatrix<T>& Output, U Scalar, const BlockMatrix<T>& Input1, const BlockMatrix<T>& Input2 )
{
    assert(Output.get_nblocks()!=0);

    Output.fill(0.);
    auto& ci = MeshGrid<R>::ConvolutionIndex[{Output.meshgrid->get_id(), Input1.meshgrid->get_id(), Input2.meshgrid->get_id()}];

    if(ci.get_Size(0) == 0 ){
        MeshGrid<R>::Calculate_ConvolutionIndex(*(Output.meshgrid), *(Input1.meshgrid), *(Input2.meshgrid));

    }
    for(int iblock_o=0; iblock_o<Output.get_nblocks(); iblock_o++){
        for(int iblock_i2=0; iblock_i2<Input2.get_nblocks(); iblock_i2++){
            auto& iblock_i1 = ci(iblock_o, iblock_i2);
            Matrix_gemm(Output[iblock_o], Scalar+im*0., Input1[iblock_i1], Input2[iblock_i2], 1.+im*0.);
        }
    }
}


template<typename T>
BlockMatrix<T>::~BlockMatrix()
{
    Values.~mdarray<T,3>();
}