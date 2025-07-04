
#ifndef VECTOR_HPP
#define VECTOR_HPP

#include <complex>
#include <memory>
#include <cassert>
#define MKL_Complex16 std::complex<double>
#include "mkl.h"
#include <iostream>
#include <iomanip>
#include <type_traits>

#include <cstdarg>

#include "mdContainers/mdContainers.hpp"
#include "LinearAlgebra/gemm.hpp"

template<typename T>
class Matrix;

template<typename T>
class Vector{
    private:
        mdarray<T,1> Values;
    public:
        Vector() = default;

        Vector(const Vector& v){this->Values = v.Values; };
        Vector& operator=(const Vector& v){this->Values = v.Values; return *this;};
        
        Vector(Vector&& v) = default;
        Vector& operator=(Vector&& v) = default;

        Vector(mdarray<T,1> Values_) : Values(std::move(Values_)){};
        
        Vector(const int& n);
        void initialize(const int& n);

        template <typename... Args>
        void initialize_n(Args... args)
        {
            int dim = sizeof...(args);
            std::vector<std::common_type_t<Args...>> aux{args...};
            Values.initialize({dim});
            for(int ix=0; ix<dim; ix++){
                Values[ix] = aux[ix];
            }
        }

        const T& operator()(const int& n) const;
        T& operator()(const int& n);
        const T& operator[](const int& n) const;
        T& operator[](const int& n);

        //Matrix<T> operator*(const Matrix<T>& B) const;
        Vector<T> operator*(T Scalar) const;
        
        template<typename U>
        friend Vector<U> operator*(U Scalar, const Vector<U> v);

        Vector operator-() const;
        Vector operator+() const;
        Vector operator-(const Vector& v) const;
        Vector operator+(const Vector& v) const;
        Vector& operator+=(const Vector& v);
        Vector& operator-=(const Vector& v);

        double norm() const;

        const T* data() const {return Values.data();};
        T* data() {return const_cast<T*>((static_cast<const Vector<T>&>(*this)).Values.data());};

        auto begin() const { return Values.begin(); }
        auto end() const {return Values.end(); } 	
        int get_NumberOfElements() const;
        

        friend class Matrix<T>;
};


//overloading writing matrix
template<class T>
std::ostream& operator<<(std::ostream& os, const Vector<T>& v);
#include "Vector_definitions.hpp"

#endif
