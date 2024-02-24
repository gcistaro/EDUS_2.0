#ifndef BLOCKMATRIX_HPP
#define BLOCKMATRIX_HPP

#include <complex>
#include <memory>
#include <cassert>
#define MKL_Complex16 std::complex<double>
#include "mkl.h"
#include <iostream>
#include <iomanip>
#include <typeinfo>

#include "mdContainers/mdContainers.hpp"
#include "MeshGrid/MeshGrid.hpp"
#include "Geometry/Matrix.hpp"



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

        std::vector<Matrix<T>> submatrix;
        std::shared_ptr<MeshGrid<space>> meshgrid;    
        Matrix<T> EmptyMatrix;    
        T NullValue = 0;
        void initialize_submatrix();
    public:
        mdarray<T,3> Values; //first index -> block, others -> matrix
        BlockMatrix() : Values(mdarray<T,3>()), EmptyMatrix(Matrix<T>()){};
        BlockMatrix(const size_t& nblocks, const size_t& nrows, const size_t& ncols){this->initialize(nblocks, nrows, ncols);}
        void initialize(mdarray<T,3>& Values);
        void initialize(const size_t& nblocks, const size_t& nrows, const size_t& ncols);
        BlockMatrix(const BlockMatrix<T, space>& A);
        BlockMatrix<T, space>& operator=(const BlockMatrix<T, space>& m);
        
        BlockMatrix(BlockMatrix<T, space>&& A);
        BlockMatrix<T, space>& operator=(BlockMatrix<T, space>&& m);

        const T& operator()(const int& nk1, const int& mk2) const; //full index
        T& operator()(const int& nk1, const int& mk2);             //full index

        const T& operator()(const int& iblock, const int& n, const int& m) const; //block-like index
        T& operator()(const int& iblock, const int& n, const int& m);             //block-like index

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
	    std::shared_ptr<MeshGrid<space>>& get_MeshGrid(){return this->meshgrid;};

        void set_MeshGrid(const MeshGrid<space>& meshgrid_){std::cout << "Gonna create shared ptr\n"; this->meshgrid = std::make_shared<MeshGrid<space>>(meshgrid_); std::cout << "DONE!\n";};
        bool has_meshgrid()
        {
            return (*meshgrid != nullptr);
        }
        
        template<typename T_, Space space_, typename U>
        friend void convolution(BlockMatrix<T_,space_>& Output, U Scalar, const BlockMatrix<T_,space_>& Input1, const BlockMatrix<T_,space_>& Input2 );
        
        void diagonalize(std::vector<mdarray<double,1>>& Eigenvalues,
                         BlockMatrix<std::complex<double>,space>& Eigenvectors) const;


        void test_submatrix();

        //destructor
        ~BlockMatrix();	
};




//overloading writing matrix
template<class T>
std::ostream& operator<<(std::ostream& os, BlockMatrix<T>& m);
#include "BlockMatrix_definitions.hpp"

#endif
