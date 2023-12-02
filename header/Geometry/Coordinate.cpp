void Basis::initialize(const Matrix<double>& A)
{
	std::cout << "M=A" << std::endl;
    M=A;
    std::cout << "invM \n";
    invM=M.inverse();
    std::cout << "MetricTensor \n";
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
Coordinate<space>::Coordinate(const double& c1, const double& c2, const double& c3)
{
    this->initialize(c1,c2,c3);
}

template<Space space>
Coordinate<space>::Coordinate(const double& c1, const double& c2, const double& c3, const std::string& key)
{
    this->initialize(c1,c2,c3,key);
}

template<Space space>
void Coordinate<space>::initialize(const double& c1, const double& c2, const double& c3)
{
    this->initialize(c1,c2,c3,std::string("Cartesian"));
}

template<Space space>
void Coordinate<space>::initialize(const double& c1, const double& c2, const double& c3, const std::string& KeyForBasis)
{
    if(KeyForBasis != "Cartesian"){
        assert(BasisDictionary.find(KeyForBasis) != BasisDictionary.end());
    }
    CoordinateDictionary[KeyForBasis] = Vector<double>();
    CoordinateDictionary[KeyForBasis].initialize_n(c1,c2,c3);
    if(KeyForBasis != std::string("Cartesian")){
        auto& M = BasisDictionary[KeyForBasis].get_invM();
        CoordinateDictionary["Cartesian"] = M*CoordinateDictionary[KeyForBasis];
    }
}

template<Space space>
Coordinate<space>::Coordinate(const Coordinate<space>& v)
{
    *this = v;
}

template<Space space>
Coordinate<space>& Coordinate<space>::operator=(const Coordinate<space>& v)
{
    this->CoordinateDictionary = v.CoordinateDictionary;
    return *this;
}

template<Space space>
Coordinate<space>::Coordinate(Coordinate<space>&& v)
{
    *this = v;
}

template<Space space>
Coordinate<space>& Coordinate<space>::operator=(Coordinate<space>&& v)
{
    this->CoordinateDictionary = std::move(v.CoordinateDictionary);
    return *this;
}

template<Space space>
Coordinate<space> Coordinate<space>::operator-() const
{
    return(Coordinate<space>(-CoordinateDictionary.at("Cartesian")[0],
                             -CoordinateDictionary.at("Cartesian")[1],
                             -CoordinateDictionary.at("Cartesian")[2]));
}

template<Space space>
Coordinate<space> Coordinate<space>::operator-(const Coordinate<space>& v) const
{
    return (*this)+(-v);
}

template<Space space>
Coordinate<space> Coordinate<space>::operator+(const Coordinate<space>& v) const
{
    Coordinate aux;
    aux.CoordinateDictionary["Cartesian"] = this->CoordinateDictionary.at("Cartesian") + 
                                            v.CoordinateDictionary.at("Cartesian"); 
    return aux;
}

template<Space space>
Coordinate<space>& Coordinate<space>::operator+=(const Coordinate& v)
{
    (*this) = (*this)+v;
    return (*this);
}


template<Space space>
Coordinate<space> Coordinate<space>::operator*(const double& alpha) const
{
    Coordinate<space> v1;
    v1.CoordinateDictionary["Cartesian"] = alpha*this->CoordinateDictionary.at("Cartesian");
    return v1;
}

template<Space space>
Coordinate<space> operator*(const double& alpha, const Coordinate<space>& v)
{
    return v*alpha;
}



template<Space space>
std::ostream& operator<<(std::ostream& os, const Coordinate<space>& v)
{
    os << v.CoordinateDictionary;
    return os;
}


template<Space space>
void Coordinate<space>::add_Coordinate(const std::string& KeyForBasis)
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

template<Space space>
const Vector<double>& Coordinate<space>::get(const std::string& KeyForBasis) const
{
    (const_cast<Coordinate<space>&>(*this)).add_Coordinate(KeyForBasis);
    return CoordinateDictionary.at(KeyForBasis);
}

template<Space space>
Vector<double>& Coordinate<space>::get(const std::string& KeyForBasis)
{
    this->add_Coordinate(KeyForBasis);
    return CoordinateDictionary[KeyForBasis];
}

template<Space space>
double Coordinate<space>::norm() const
{
    auto& cart = CoordinateDictionary.at("Cartesian");
    return (std::sqrt(cart[0]*cart[0]+cart[1]*cart[1]+cart[2]*cart[2]));
}

