#ifndef GEMM_HPP
#define GEMM_HPP


#include <type_traits>
#include "mkl.h"

template<typename U, typename V>
void gemm(int m, int n1, int k1, U alpha, const double* A, int k2, 
                   const double* B, int n2, V beta, double* C, int n3)
{
    //if(std::is_same<T,double>::value){
    //    double alpha_ = alpha;
    //    double beta_ = beta;
        
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
                    m, n1, k1, alpha, A, k2, B, n2, beta, C, n3);

    //}
}

template<typename U, typename V>
void gemm(int m, int n1, int k1, U alpha, const std::complex<double>* A, int k2, 
                   const std::complex<double>* B, int n2, V beta, std::complex<double>* C, int n3)
{
    //else if(std::is_same<T,std::complex<double>>::value){
    void* alpha_ = &alpha;
    void* beta_ = &beta;

        cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
                    m, n1, k1, alpha_, A, k2, B, n2, beta_, C, n3);
    //}
}

#endif