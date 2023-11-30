#ifndef BLOCKMATRIX_HPP
#define BLOCKMATRIX_HPP

#include <complex>
#include <memory>
#include <cassert>
#include "mkl.h"
#include <iostream>
#include <iomanip>
#include <typeinfo>

#include "../mdContainers/mdContainers.hpp"
#include "../MeshGrid/MeshGrid.hpp"
#include "../Geometry/Matrix.hpp"



//Our sparse matrices are block matrices. 
//we save them in a contiguous array of size (#blocks, #rows, #cols)
//if the blocks have different shape, #rows and #cols must be the maximum. 
//in the current rudimental implementation, the blocks are all the same and 
//to map matrix elements (nk, n'k') to indices of the 3d array Values
//we have an auxiliary vector. The matrix element is non zero only if
//its 1d index is between SubMatrixPointer.first and SubMatrixPointer.first+nrows*ncols
//for one SubMatrixPointer. 
//we need two different styles
//1 - O(nk, nk') = O(n,n')(k) delta(k-k')
//2 - O(nR, nR') = O(n,n')(R-R') 
template<typename T=std::complex<double>, Space space=R>
class BlockMatrix{
    private:
        mdarray<T,3> Values; //first index -> block, others -> matrix
        std::vector<Matrix<T>> submatrix;
        std::shared_ptr<MeshGrid<space>> meshgrid;    
        Matrix<T> EmptyMatrix;    
    public:
        BlockMatrix(){};
        BlockMatrix(const int& nblocks, const int& nrows, const int& ncols);
        void initialize(mdarray<T,3>& Values);

        BlockMatrix(const BlockMatrix<T>& A);
        BlockMatrix& operator=(const BlockMatrix<T>& m);
        
        BlockMatrix(BlockMatrix<T>&& A);
        BlockMatrix<T>& operator=(BlockMatrix<T>&& m);

        //const T& operator()(const int& nk1, const int& mk2) const; //full index
        //T& operator()(const int& nk1, const int& mk2);             //full index

        //const T& operator()(const int& iblock, const int& n, const int& m) const; //block-like index
        //T& operator()(const int& iblock, const int& n, const int& m);             //block-like index

        void fill(const T& Scalar);

        Matrix<T>& operator[](const int& iblock);
        const Matrix<T>& operator[](const int& iblock) const;

        BlockMatrix<T> operator*(const BlockMatrix<T>& B);

        const T* data() const {return Values.data();};
        T* data() {return const_cast<T*>((static_cast<const BlockMatrix<T,space>&>(*this)).Values.data());};
        //friend void multiply(Matrix<T>& OutputMatrix, const auto& Scalar1, const Matrix<T>& Matrix1, 
        //                                              const auto& Scalar2, const Matrix<T>& Matrix2);

        auto begin() const { return Values.begin(); }
        auto end() const {return Values.end(); } 	
        size_t get_nblocks() const{return Values.get_Size(0);};
        size_t get_nrows() const{return Values.get_Size(1);};
        size_t get_ncols() const{return Values.get_Size(2);};
        MeshGrid<space>* get_MeshGrid()const {return this->meshgrid.get();};

        void set_MeshGrid(MeshGrid<space>& meshgrid_){ this->meshgrid = std::make_shared<MeshGrid<space>>(meshgrid_);};
        template<typename T_, Space space_, typename U>
        friend void convolution(BlockMatrix<T_,space_>& Output, U Scalar, const BlockMatrix<T_,space_>& Input1, const BlockMatrix<T_,space_>& Input2 );
        
	//destructor
        ~BlockMatrix();	
};




//overloading writing matrix
template<class T>
std::ostream& operator<<(std::ostream& os, BlockMatrix<T>& m);
#include "BlockMatrix.cpp"

#endif
