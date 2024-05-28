#include <string>
#include <map>
#include "Constants.hpp"
#include "Vector.hpp"
#include "Matrix.hpp"

#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

template<typename T>
using Dictionary = std::map<std::string, T>;//WARNING !! operator[] works only if the Dictionary is not const.


template<typename T>
std::ostream& operator<<(std::ostream& os, const Dictionary<T>& dic)
{
    for(const auto& [key,value]: dic){
        os << "key: " << key << std::endl;
        os << value;
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
        
        inline const Matrix<double>& get_M() const {return M;};
        inline const Matrix<double>& get_invM() const {return invM;};
        inline const Matrix<double>& get_MetricTensor() const {return MetricTensor;};

        friend std::ostream& operator<<(std::ostream& os, const Basis& basis);
};

std::ostream& operator<<(std::ostream& os, const Basis& basis);



class Coordinate{
    private:
        Dictionary<Vector<double>> CoordinateDictionary;
        static Dictionary<Basis> BasisDictionary;

    public:
        Coordinate(){};
        Coordinate(const double& c1, const double& c2, const double& c3);
        Coordinate(const double& c1, const double& c2, const double& c3, const std::string& KeyForBasis);
        Coordinate(const Coordinate& v);
        Coordinate& operator=(const Coordinate& v);

        Coordinate(Coordinate&& v);
        Coordinate& operator=(Coordinate&& v);

        void initialize(const double& c1, const double& c2, const double& c3);
        void initialize(const double& c1, const double& c2, const double& c3, const std::string& KeyForBasis);
        inline static void add_Basis(const Basis& Basis_to_add, const std::string& KeyForBasis);
        inline static Basis& get_Basis(const std::string& KeyForBasis){
            return BasisDictionary[KeyForBasis];
        }

        Coordinate operator-() const;
        Coordinate operator-(const Coordinate& v) const;
        Coordinate operator+(const Coordinate& v) const;
        Coordinate& operator+=(const Coordinate& v);

        Coordinate operator*(const double& alpha) const;
        Coordinate operator/(const double& alpha) const;

        friend Coordinate operator*(const double& alpha, const Coordinate& v);

        void add_Coordinate(const std::string& KeyForBasis);
        inline const Vector<double>& get(const std::string& KeyForBasis) const;
        inline Vector<double>& get(const std::string& KeyForBasis);
        inline double norm() const;
        inline double dot(const Coordinate& v1) const;
        
        friend std::ostream& operator<<(std::ostream& os, const Coordinate& v);

};  


//inline functions defined here in hpp
void Coordinate::add_Basis(const Basis& Basis_to_add, const std::string& KeyForBasis)
{
    BasisDictionary[KeyForBasis] = Basis_to_add;
}

const Vector<double>& Coordinate::get(const std::string& KeyForBasis) const
{
    (const_cast<Coordinate&>(*this)).add_Coordinate(KeyForBasis);
    return CoordinateDictionary.at(KeyForBasis);
}

Vector<double>& Coordinate::get(const std::string& KeyForBasis)
{
    this->add_Coordinate(KeyForBasis);
    return CoordinateDictionary[KeyForBasis];
}

double Coordinate::norm() const
{
    auto& cart = CoordinateDictionary.at("Cartesian");
    return (std::sqrt(cart(0)*cart(0)+cart(1)*cart(1)+cart(2)*cart(2)));
}

double Coordinate::dot(const Coordinate& v) const 
{
    auto& x = this->get("Cartesian");
    auto& y = v.get("Cartesian");
    return ( x[0] * y[0] + x[1] * y[1] + x[2] * y[2] );
}

#endif
