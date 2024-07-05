template<class T>
Vector<T>::Vector(const size_t& n)
{
    this->initialize(n);
}

template<class T>
void Vector<T>::initialize(const size_t& n)
{
    (*this).~Vector<T>();
    this->Values.initialize({n});
}


//copy constructor
//template<class T>
//Vector<T>::Vector(const Vector<T>& v){
//    *this = v;
//}
//
//
////copy assigment
//template<class T>
//Vector<T>& Vector<T>::operator=(const Vector<T>& v){
//    (*this).~Vector<T>();
//    this->Values = v.Values;
//    return *this;
//}
//
//
////move constructor
//template<class T>
//Vector<T>::Vector(Vector&& v){
//    *this = v;
//}
//
//
////move assignment
//template<class T>
//Vector<T>& Vector<T>::operator=(Vector<T>&& v){
//    (*this).~Vector<T>();
//    (*this).Values = v.Values;
//    return *this;
//}


template<class T>
const T& Vector<T>::operator()(const int& i) const
{
    assert ( i >= 0 && i < this->get_NumberOfElements() );
    return this->Values(i);
}


template<class T>
T& Vector<T>::operator()(const int& i) 
{
    return const_cast<T&>(static_cast<const Vector<T>&> (*this)(i));
}

template<class T>
const T& Vector<T>::operator[](const int& i) const
{
    return this->Values(i);
}


template<class T>
T& Vector<T>::operator[](const int& i) 
{
    return this->Values(i);
}


template<class T>
Vector<T> Vector<T>::operator*(T Scalar) const
{
    Vector<T> C(this->get_NumberOfElements());
    for(int i=0; i<this->get_NumberOfElements(); i++){
        C(i) = Scalar*(*this)(i);
    }
    return C;
}

template<class T>
Vector<T> operator*(T Scalar, const Vector<T> A)
{
    return A*Scalar;
}

template<class T>
Vector<T> Vector<T>::operator-() const
{
    Vector<T> v1(this->get_NumberOfElements());
    for(int ix=0; ix<v1.get_NumberOfElements(); ix++){
        v1(ix) = -(*this)(ix);
    }
    return v1;
}

template<class T>
Vector<T> Vector<T>::operator+() const
{
    return *this;
}

template<class T>
Vector<T> Vector<T>::operator+(const Vector& v) const
{
    assert(this->get_NumberOfElements() == v.get_NumberOfElements());
    Vector<T> v1(this->get_NumberOfElements());
    
    for(int ix=0; ix<v1.get_NumberOfElements(); ix++){
        v1(ix) = (*this)(ix) + v(ix);    
    }
    return v1;
}

template<class T>
Vector<T> Vector<T>::operator-(const Vector& v) const
{
    return (*this) + (-v);
}

template<class T>
Vector<T>& Vector<T>::operator+=(const Vector& v)
{
    (*this) = (*this)+v;
    return (*this);
}

template<class T>
Vector<T>& Vector<T>::operator-=(const Vector& v)
{
    return (*this) += (-v);
}


template<class T>
size_t Vector<T>::get_NumberOfElements() const
{
    return this->Values.get_Size(0);
}


//destructor
template<typename T>
Vector<T>::~Vector()
{
    (*this).Values.~mdarray<T,1>();
}

template<class T>
double Vector<T>::norm() const
{
    double norm = 0.;
    for( auto& v : *this ) {
        norm += std::abs(v)*std::abs(v);
    }
    return std::sqrt(norm);
}


template<class T>
std::ostream& operator<<(std::ostream& os, const Vector<T>& v)
{
    for(int i=0; i<v.get_NumberOfElements(); i++){
        os << std::setprecision(15) << std::setw(24) << v(i); 
    }
    os << std::endl;
    return os;
}
