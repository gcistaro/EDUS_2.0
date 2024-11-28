#include <complex>
#include "Operator/Operator.hpp"
#include "DESolver.hpp"


template<>
void DESolver<Operator<std::complex<double>>>::initialize(Operator<std::complex<double>>& Function_, 
                const std::function<void(Operator<std::complex<double>>&)>& EvaluateInitialCondition_, 
                const std::function<void(Operator<std::complex<double>>&, 
                const double&, const Operator<std::complex<double>>&)>& EvaluateSourceFunction_,
                SolverType type__)
{
    type = type__;
    Function = &Function_;
    EvaluateInitialCondition = EvaluateInitialCondition_;
    EvaluateSourceFunction = EvaluateSourceFunction_;
    EvaluateInitialCondition(*Function);
    //need to initialize fft!
    for(auto& aux_f : aux_Function) {
        aux_f.initialize_fft(*(Function));
        aux_f.lock_space(Space::k);
    }
}