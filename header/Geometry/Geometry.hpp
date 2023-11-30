#include <string>
#include <map>
#include "Matrix.hpp"

#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

template<typename T>
using Dictionary = std::map<std::string, T>;//WARNING !! operator[] works only if the Dictionary is not const.

enum Space{k,R};

template<typename T>
std::ostream& operator<<(std::ostream& os, const Dictionary<T>& dic)
{
    for(const auto& [key,value]: dic){
        std::cout << "key: " << key << std::endl;
        std::cout << value;
    }
    return os;
}


class Basis{
    private:
        Matrix<double> M;
        Matrix<double> invM;
        Matrix<double> MetricTensor;
        Matrix<double> InverseMetricTensor;
    public:
        Basis(){M.initialize(3,3); invM.initialize(3,3); MetricTensor.initialize(3,3); InverseMetricTensor.initialize(3,3);};
        Basis(const Matrix<double>& A) {this->initialize(A);};
        void initialize(const Matrix<double>& A);
        
        inline const Matrix<double>& get_M() const;
        inline const Matrix<double>& get_invM() const;
        inline const Matrix<double>& get_MetricTensor() const;

        friend std::ostream& operator<<(std::ostream& os, const Basis& basis);
};

std::ostream& operator<<(std::ostream& os, const Basis& basis)
{
    os<< "Basis"; 
    os<< "\n        M:\n";
    os<< basis.M;
    os<< "\n     invM:\n";
    os<< basis.invM;

    return os;
}


class Coordinate
{
    private:
        std::array<double,3> Values{-1.e+08,-1.e+08, -1.e08};
    public:
        Coordinate(){};
        Coordinate(const double& c1, const double& c2, const double& c3) : Values({c1,c2,c3}){};

        double& operator[](const int& i){return const_cast<double&>(static_cast<const Coordinate&>(*this)[i]);};
        const double& operator[](const int& i) const {return Values[i];};
        
        template<Space space> friend class Vector;

        friend std::ostream& operator<<(std::ostream& os, const Coordinate& c_);
};

std::ostream& operator<<(std::ostream& os, const Coordinate& c_)
{
    std::cout << c_.Values[0] << " " << c_.Values[1] << " " << c_.Values[2] << "\n";
    return os;
}


template<Space space>
class Vector{
    private:
        Dictionary<Coordinate> CoordinateDictionary;
        static Dictionary<Basis> BasisDictionary;

    public:
        Vector(){};
        Vector(const double& c1, const double& c2, const double& c3);
        Vector(const double& c1, const double& c2, const double& c3, const std::string& key);
        Vector(const Vector& v);
        Vector& operator=(const Vector& v);

        Vector(Vector&& v);
        Vector& operator=(Vector&& v);

        void initialize(const double& c1, const double& c2, const double& c3);
        void initialize(const double& c1, const double& c2, const double& c3, const std::string& key);
        inline static void add_Basis(const Basis& Basis_to_add, const std::string& key);
        inline static Basis& get_Basis(const std::string& KeyForBasis){
            return BasisDictionary[KeyForBasis];
        }

        Vector operator-() const;
        Vector operator-(const Vector& v) const;
        Vector operator+(const Vector& v) const;
        Vector& operator+=(const Vector& v);

        Vector operator*(const double& alpha) const;

        template<Space space_>
        friend Vector<space_> operator*(const double& alpha, const Vector<space_>& v);

        void add_Coordinate(const std::string& KeyForBasis);
        inline const Coordinate& get(const std::string& KeyForBasis) const;
        inline Coordinate& get(const std::string& KeyForBasis);
        inline double norm() const;
        
        template<Space space_>
        friend std::ostream& operator<<(std::ostream& os, const Vector<space_>& v);

};  


template<Space space>
Dictionary<Basis> Vector<space>::BasisDictionary;


template<Space space>
void Vector<space>::add_Basis(const Basis& Basis_to_add, const std::string& key)
{
    BasisDictionary[key] = Basis_to_add;
}



#include "Geometry.cpp"


#endif
