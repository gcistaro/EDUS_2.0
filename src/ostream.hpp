#include <iostream>
#include <complex>
#include <iomanip>
#include <vector>

std::ostream& operator<<(std::ostream& os, const std::complex<double>& number);
std::ostream& operator<<(std::ostream& os, const std::vector<std::complex<double>>& vector);