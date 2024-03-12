#include "Coordinate.hpp"
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
