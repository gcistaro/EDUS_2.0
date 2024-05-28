
#include "ostream.hpp"

std::ostream& operator<<(std::ostream& os, const std::complex<double>& number)
{
    os << std::setw(15) << std::setprecision(8) << number.real();
    os << std::setw(15) << std::setprecision(8) << number.imag();
    return os;
}


std::ostream& operator<<(std::ostream& os, const std::vector<std::complex<double>>& vector)
{
    for(auto& number : vector) {
        os << number;
    }
    return os;
}