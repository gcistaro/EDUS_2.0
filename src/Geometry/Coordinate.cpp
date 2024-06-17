#include "Coordinate.hpp"
Dictionary<Basis> Coordinate::BasisDictionary;

void Basis::initialize(const Matrix<double>& A)
{
    M=A;
    invM=M.inverse();
    MetricTensor = M.transpose()*M;
    InverseMetricTensor = invM.transpose()*invM;
}

std::ostream& operator<<(std::ostream& os, const Basis& basis)
{
    os<< "Basis"; 
    os<< "\n        M:\n";
    os<< basis.M;
    os<< "\n     invM:\n";
    os<< basis.invM;

    return os;
}

Coordinate::Coordinate(const double& c1, const double& c2, const double& c3)
{
    this->initialize(c1,c2,c3);
}

Coordinate::Coordinate(const double& c1, const double& c2, const double& c3, const std::string& key)
{
    this->initialize(c1,c2,c3,key);
}

void Coordinate::initialize(const double& c1, const double& c2, const double& c3)
{
    this->initialize(c1,c2,c3,std::string("Cartesian"));
}

void Coordinate::initialize(const double& c1, const double& c2, const double& c3, const std::string& KeyForBasis)
{
    if(KeyForBasis != "Cartesian"){
        assert(BasisDictionary.find(KeyForBasis) != BasisDictionary.end());
    }
    CoordinateDictionary[KeyForBasis] = Vector<double>();
    CoordinateDictionary.at(KeyForBasis).initialize_n(c1,c2,c3);
    if(KeyForBasis != std::string("Cartesian")){
        auto& M = BasisDictionary[KeyForBasis].get_M();
        CoordinateDictionary["Cartesian"] = M*CoordinateDictionary[KeyForBasis];
    }
}

Coordinate::Coordinate(const Coordinate& v)
{
    *this = v;
}

Coordinate& Coordinate::operator=(const Coordinate& v)
{
    this->CoordinateDictionary = v.CoordinateDictionary;
    return *this;
}

Coordinate::Coordinate(Coordinate&& v)
{
    *this = v;
}

Coordinate& Coordinate::operator=(Coordinate&& v)
{
    this->CoordinateDictionary = std::move(v.CoordinateDictionary);
    return *this;
}

Coordinate Coordinate::operator-() const
{
    auto& cart = CoordinateDictionary.at("Cartesian");
    return(Coordinate(-cart(0),
                             -cart(1),
                             -cart(2)));
}

Coordinate Coordinate::operator-(const Coordinate& v) const
{
    return (*this)+(-v);
}

Coordinate Coordinate::operator+(const Coordinate& v) const
{
    Coordinate aux;
    aux.CoordinateDictionary["Cartesian"] = this->CoordinateDictionary.at("Cartesian") + 
                                            v.CoordinateDictionary.at("Cartesian"); 
    return aux;
}

Coordinate& Coordinate::operator+=(const Coordinate& v)
{
    (*this) = (*this)+v;
    return (*this);
}

Coordinate Coordinate::operator*(const double& alpha) const
{
    Coordinate v1;
    v1.CoordinateDictionary["Cartesian"] = alpha*this->CoordinateDictionary.at("Cartesian");
    return v1;
}

Coordinate Coordinate::operator/(const double& alpha) const
{
    Coordinate v1;
    v1.CoordinateDictionary["Cartesian"] = this->CoordinateDictionary.at("Cartesian");

    v1.CoordinateDictionary["Cartesian"][0] = v1.CoordinateDictionary["Cartesian"][0]/alpha; 
    v1.CoordinateDictionary["Cartesian"][1] = v1.CoordinateDictionary["Cartesian"][1]/alpha; 
    v1.CoordinateDictionary["Cartesian"][2] = v1.CoordinateDictionary["Cartesian"][2]/alpha; 
    return v1;
}

Coordinate operator*(const double& alpha, const Coordinate& v)
{
    return v*alpha;
}

std::ostream& operator<<(std::ostream& os, const Coordinate& v)
{
    os << v.CoordinateDictionary;
    return os;
}

void Coordinate::add_Coordinate(const std::string& KeyForBasis)
{
    if(KeyForBasis != "Cartesian"){
        assert(BasisDictionary.find(KeyForBasis) != BasisDictionary.end());
    }

    if(CoordinateDictionary.find(KeyForBasis) == CoordinateDictionary.end())
    {
        auto& cart = CoordinateDictionary.at("Cartesian");
        auto& invM = BasisDictionary[KeyForBasis].get_invM();
        CoordinateDictionary[KeyForBasis] = invM*cart;
    }

}

