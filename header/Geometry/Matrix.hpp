#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <complex>
#include <memory>
#include <cassert>
#include "mkl.h"
#include <iostream>
#include <iomanip>
#include <type_traits>
#include <cassert>
#include "../mdContainers/mdContainers.hpp"
#include "../LinearAlgebra/gemm.hpp"

template<typename T>
class Vector;


template<class T>
class Matrix{
    private:
        mdarray<T,2> Values;
    public:
        Matrix(){Values = mdarray<T,2>();};
        Matrix(mdarray<T,2> Values_) : Values(std::move(Values_)){};
        Matrix(const int& nrows, const int& ncols);
        void initialize(const int& nrows, const int& ncols);

        Matrix(const Matrix& A);
        Matrix& operator=(const Matrix& m);
        
        Matrix(Matrix&& A);
        Matrix& operator=(Matrix&& m);

        const T& operator()(const int& n, const int& m) const;
        T& operator()(const int& n, const int& m);

        Matrix<T> operator*(const Matrix<T>& B) const;
        Matrix<T> operator*(T Scalar) const;
        Vector<T> operator*(const Vector<T>& v) const;
         
        template<typename U>
        friend Matrix<U> operator*(U Scalar, const Matrix<U> M);

        Matrix<T> inverse() const;
        Matrix<T> transpose() const;
        void LUdecompose(Matrix<T>& LU, lapack_int** pointer_to_ipiv) const;
        //Matrix<T> LUdecompose() const;
	T determinant() const;
        const T* data() const {return Values.data();};
        T* data() {return const_cast<T*>((static_cast<const Matrix<T>&>(*this)).Values.data());};
        //friend void multiply(Matrix<T>& OutputMatrix, const auto& Scalar1, const Matrix<T>& Matrix1, 
        //                                              const auto& Scalar2, const Matrix<T>& Matrix2);

        auto begin() const { return Values.begin(); }
        auto end() const {return Values.end(); } 	
        int get_nrows() const;
        int get_ncols() const;
        
	//destructor
        ~Matrix();	

        friend class Vector<T>;
};


template<typename T>
void Matrix_gemm(Matrix<T>& OutputMatrix, const T& alpha, const Matrix<T>& InputMatrix1, const Matrix<T>& InputMatrix2, const T& beta);

//overloading writing matrix
template<class T>
std::ostream& operator<<(std::ostream& os, const Matrix<T>& m);
#include "Matrix.cpp"

#endif
