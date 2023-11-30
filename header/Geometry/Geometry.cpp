void Basis::initialize(const Matrix<double>& A)
{
    M=A;
    invM=M.inverse();
    MetricTensor = M.transpose()*M;
    InverseMetricTensor = invM.transpose()*invM;
}



const Matrix<double>& Basis::get_M() const
{
    return M;
}

const Matrix<double>& Basis::get_invM() const 
{
    return invM;
}

const Matrix<double>& Basis::get_MetricTensor() const 
{
    return MetricTensor;
}

template<Space space>
Vector<space>::Vector(const double& c1, const double& c2, const double& c3)
{
    this->initialize(c1,c2,c3);
}

template<Space space>
Vector<space>::Vector(const double& c1, const double& c2, const double& c3, const std::string& key)
{
    this->initialize(c1,c2,c3,key);
}

template<Space space>
void Vector<space>::initialize(const double& c1, const double& c2, const double& c3)
{
    this->initialize(c1,c2,c3,std::string("Cartesian"));
}

template<Space space>
void Vector<space>::initialize(const double& c1, const double& c2, const double& c3, const std::string& KeyForBasis)
{
    if(KeyForBasis != "Cartesian"){
        assert(BasisDictionary.find(KeyForBasis) != BasisDictionary.end());
    }
    CoordinateDictionary[KeyForBasis] = Coordinate(c1,c2,c3);
    if(KeyForBasis != std::string("Cartesian")){
        std::array<double,3> cart{0};
        //following line means: cart = BasisDictionary[key].M*CoordinateDictionary[key]
        gemm(3,1,3,1.,&(BasisDictionary[KeyForBasis].get_M()(0,0)),3,
             &(CoordinateDictionary[KeyForBasis][0]),1,0.,&cart[0],1);
        CoordinateDictionary[std::string("Cartesian")] = Coordinate({cart[0],cart[1],cart[2]});
    }
}

template<Space space>
Vector<space>::Vector(const Vector<space>& v)
{
    *this = v;
}

template<Space space>
Vector<space>& Vector<space>::operator=(const Vector<space>& v)
{
    this->CoordinateDictionary = v.CoordinateDictionary;
    return *this;
}

template<Space space>
Vector<space>::Vector(Vector<space>&& v)
{
    *this = v;
}

template<Space space>
Vector<space>& Vector<space>::operator=(Vector<space>&& v)
{
    this->CoordinateDictionary = std::move(v.CoordinateDictionary);
    return *this;
}

template<Space space>
Vector<space> Vector<space>::operator-() const
{
    return(Vector<space>(-CoordinateDictionary.at("Cartesian")[0], 
                         -CoordinateDictionary.at("Cartesian")[1], 
                         -CoordinateDictionary.at("Cartesian")[2]));
}

template<Space space>
Vector<space> Vector<space>::operator-(const Vector<space>& v) const
{
    return (*this)+(-v);
}

template<Space space>
Vector<space> Vector<space>::operator+(const Vector<space>& v) const
{
    Vector aux;
    for(int ix=0; ix<3; ix++){
        aux.CoordinateDictionary[std::string("Cartesian")][ix] = 
            this->CoordinateDictionary.at("Cartesian")[ix] + v.CoordinateDictionary.at("Cartesian")[ix]; 
    }
    return aux;
}

template<Space space>
Vector<space>& Vector<space>::operator+=(const Vector& v)
{
    (*this) = (*this)+v;
    return (*this);
}


template<Space space>
Vector<space> Vector<space>::operator*(const double& alpha) const
{
    Vector<space> v1;
    for(auto icoor=0; icoor<3; icoor++){
        v1.CoordinateDictionary["Cartesian"][icoor] = alpha*this->CoordinateDictionary.at("Cartesian")[icoor];
    }
    return v1;
}

template<Space space>
Vector<space> operator*(const double& alpha, const Vector<space>& v)
{
    return v*alpha;
}



template<Space space>
std::ostream& operator<<(std::ostream& os, const Vector<space>& v)
{
    os << v.CoordinateDictionary;
    return os;
}


template<Space space>
void Vector<space>::add_Coordinate(const std::string& KeyForBasis)
{
    if(KeyForBasis != "Cartesian"){
        assert(BasisDictionary.find(KeyForBasis) != BasisDictionary.end());
    }

    if(CoordinateDictionary.find(KeyForBasis) == CoordinateDictionary.end())
    {
        auto& cart = CoordinateDictionary.at("Cartesian");
        //following line means: CoordinateDictionary[KeyForBasis] = BasisDictionary[KeyForBasis].invM*CoordinateDictionary[std::string("Cartesian")]
        gemm(3,1,3,1.,&(BasisDictionary[KeyForBasis].get_invM()(0,0)),3,
             &cart[0],1,0.,&(CoordinateDictionary[KeyForBasis][0]),1);
    }

}

template<Space space>
const Coordinate& Vector<space>::get(const std::string& KeyForBasis) const
{
    (const_cast<Vector<space>&>(*this)).add_Coordinate(KeyForBasis);
    return CoordinateDictionary.at(KeyForBasis);
}

template<Space space>
Coordinate& Vector<space>::get(const std::string& KeyForBasis)
{
    this->add_Coordinate(KeyForBasis);
    return CoordinateDictionary[KeyForBasis];
}

template<Space space>
double Vector<space>::norm() const
{
    auto& cart = CoordinateDictionary.at("Cartesian");
    return (std::sqrt(cart[0]*cart[0]+cart[1]*cart[1]+cart[2]*cart[2]));
}

