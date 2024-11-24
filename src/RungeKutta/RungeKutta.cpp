#include <complex>
#include "Operator/Operator.hpp"
#include "RungeKutta.hpp"

template<>
void RungeKutta<Operator<std::complex<double>>>::initialize(Operator<std::complex<double>>& Function_, 
                const std::function<void(Operator<std::complex<double>>&)>& EvaluateInitialCondition_, 
                const std::function<void(Operator<std::complex<double>>&, 
                const double&, const Operator<std::complex<double>>&)>& EvaluateSourceFunction_)
{
    Function = &Function_;
    EvaluateInitialCondition = EvaluateInitialCondition_;
    EvaluateSourceFunction = EvaluateSourceFunction_;
    EvaluateInitialCondition(*Function);
    //need to initialize fft!
    AuxiliaryFunction.initialize_fft(*(Function), "AuxFunction");
    ReducingFunction.initialize_fft(*(Function), "RedFunction");
    k.initialize_fft(*(Function), "k");
    AuxiliaryFunction.lock_space(Space::k);
    ReducingFunction.lock_space(Space::k);
    k.lock_space(Space::k);
}
