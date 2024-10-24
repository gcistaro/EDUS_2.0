#ifndef RUNGEKUTTA_HPP
#define RUNGEKUTTA_HPP
#include "DESolver/DESolver.hpp"


                                                        // -- RUNGE KUTTA DERIVED CLASS -- //

template<typename T>
class RungeKutta: public DESolver<T>{
    private:           
        T AuxiliaryFunction;   
        T ReducingFunction;    
        T k;

        using DESolver<T>::Function;
        using DESolver<T>::ResolutionTime;
        using DESolver<T>::CurrentTime;
        using DESolver<T>::EvaluateInitialCondition;
        using DESolver<T>::EvaluateSourceFunction;

    public:

        RungeKutta(T& Function_, const std::function<void(T&)>& EvaluateInitialCondition_, const std::function<void(T&, const double&, const T&)>& EvaluateSourceFunction_); 

        using DESolver<T>::DESolver;

        void initialize(T& Function_, const std::function<void(T&)>& EvaluateInitialCondition_, const std::function<void(T&, const double&, const T&)>& EvaluateSourceFunction_); 
        void Propagate();
};

template<typename T>
RungeKutta<T>::RungeKutta(T& Function_, const std::function<void(T&)>& EvaluateInitialCondition_, 
                        const std::function<void(T&, const double&, const T&)>& EvaluateSourceFunction_)
{
    initialize(Function_, EvaluateInitialCondition_, EvaluateSourceFunction_);
}

template<typename T>
void RungeKutta<T>::initialize(T& Function_, const std::function<void(T&)>& EvaluateInitialCondition_, const std::function<void(T&, const double&, const T&)>& EvaluateSourceFunction_)
{
    Function = &Function_;
    EvaluateInitialCondition = EvaluateInitialCondition_;
    EvaluateSourceFunction = EvaluateSourceFunction_;
    EvaluateInitialCondition(*Function);
    AuxiliaryFunction = *Function;
    ReducingFunction = *Function;
    k = *Function; 
}

template<typename T>
void RungeKutta<T>::Propagate()
{
    //k1=f(tn,yn)
    EvaluateSourceFunction(k, CurrentTime, *Function);

    //k2=f(tn+h/2,yn+h/2*k1)
    SumWithProduct(AuxiliaryFunction, 1., *Function, ResolutionTime/2., k);     // second argument of f(tn,yn) (in this case it would be yn+h/2*k1)
    SumWithProduct(ReducingFunction, 1., *Function, ResolutionTime/6., k);      // sum of the ks is reducing function
    EvaluateSourceFunction(k, CurrentTime+ResolutionTime/2., AuxiliaryFunction);

    //k3=f(tn+h/2, yn+h/2*k2)
    SumWithProduct(AuxiliaryFunction, 1., *Function, ResolutionTime/2., k);
    SumWithProduct(ReducingFunction, 1., ReducingFunction, ResolutionTime/3., k);   
    EvaluateSourceFunction(k,CurrentTime+ResolutionTime/2., AuxiliaryFunction); 

    //k4=f(tn+h,yn+h*k3)
    SumWithProduct(AuxiliaryFunction, 1., *Function, ResolutionTime, k);
    SumWithProduct(ReducingFunction, 1., ReducingFunction, ResolutionTime/3., k);  
    EvaluateSourceFunction(k, CurrentTime+ResolutionTime, AuxiliaryFunction); 

    //compute final function
    SumWithProduct(*Function, 1., ReducingFunction, ResolutionTime/6., k);  
    CurrentTime += ResolutionTime;
}



#endif