#include <complex>
#include "Operator/Operator.hpp"
#include "DESolver.hpp"

/// @brief Specialized initialization for Operator<std::complex<double>>. We initialize fft for aux_function so we allocate all the memory needed.
/// @param Function__ What will be propagated by DEsolver 
/// @param EvaluateInitialCondition__ Standard function containing the definition of the Function at InitialTime
/// @param EvaluateSourceFunction__  Standard function defining the derivative of Function, for the time propagation.
/// @param type__ Can be RK (Runge-Kutta) or AB (Adams-Bashforth)
/// @param order__ Order of the method used
template<>
void DESolver<Operator<std::complex<double>>>::initialize(Operator<std::complex<double>>& Function_, 
                const std::function<void(Operator<std::complex<double>>&)>& EvaluateInitialCondition_, 
                const std::function<void(Operator<std::complex<double>>&, 
                const double&, const Operator<std::complex<double>>&)>& EvaluateSourceFunction_,
                SolverType type__,
                int order__)
{
    order = order__;
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

    if( type == AB ) {
        initialize_beta();
    } 
}