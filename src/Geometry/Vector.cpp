#include "Vector.hpp"

template<>
std::complex<double> Vector<std::complex<double>>::dot( const Vector<std::complex<double>>& v__ ) const
{
    auto& x = *this;
    std::complex<double> dotproduct = 0;
    for ( int ix = 0; ix < this->Values.get_TotalSize(); ++ix ) {
        dotproduct += std::conj( x[ix] ) * v__[ix];
    }
    return dotproduct;
}