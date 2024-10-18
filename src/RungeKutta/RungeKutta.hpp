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

        /*
        using DESolver<T>::Function;
        using DESolver<T>::ResolutionTime;
        using DESolver<T>::CurrentTime;
        using DESolver<T>::EvaluateInitialCondition;
        using DESolver<T>::EvaluateSourceFunction;
        */

    public:
        using DESolver<T>::DESolver;
        //RungeKutta(T& Function_, const std::function<void(T&)>& EvaluateInitialCondition_, const std::function<void(T&, const double&, const T&)>& EvaluateSourceFunction_) 
        //: DESolver<T>(Function_, EvaluateInitialCondition_, EvaluateSourceFunction_){};

        void initialize(T& Function_, const std::function<void(T&)>& EvaluateInitialCondition_, const std::function<void(T&, const double&, const T&)>& EvaluateSourceFunction_) override; 
        void Propagate(); 
};

template<typename T>
void RungeKutta<T>::initialize(T& Function_, const std::function<void(T&)>& EvaluateInitialCondition_, const std::function<void(T&, const double&, const T&)>& EvaluateSourceFunction_)
{
    this->Function = &Function_;
    this->EvaluateInitialCondition = EvaluateInitialCondition_;
    this->EvaluateSourceFunction = EvaluateSourceFunction_;
    this->EvaluateInitialCondition(*this->Function);
    AuxiliaryFunction = *this->Function;
    ReducingFunction = *this->Function;
    k = *this->Function;  
}

template<typename T>
void RungeKutta<T>::Propagate()
{
    //k1=f(tn,yn)
    this->EvaluateSourceFunction(k, this->CurrentTime, *this->Function);

    //k2=f(tn+h/2,yn+h/2*k1)
    SumWithProduct(this->AuxiliaryFunction, 1., *this->Function, this->ResolutionTime/2., k);     // second argument of f(tn,yn) (in this case it would be yn+h/2*k1)
    SumWithProduct(this->ReducingFunction, 1., *this->Function, this->ResolutionTime/6., k);      // sum of the ks is reducing function
    this->EvaluateSourceFunction(k, this->CurrentTime+this->ResolutionTime/2., this->AuxiliaryFunction);

    //k3=f(tn+h/2, yn+h/2*k2)
    SumWithProduct(this->AuxiliaryFunction, 1., *this->Function, this->ResolutionTime/2., k);
    SumWithProduct(this->ReducingFunction, 1., this->ReducingFunction, this->ResolutionTime/3., k);   
    this->EvaluateSourceFunction(k,this->CurrentTime+this->ResolutionTime/2., this->AuxiliaryFunction); 

    //k4=f(tn+h,yn+h*k3)
    SumWithProduct(this->AuxiliaryFunction, 1., *this->Function, this->ResolutionTime, k);
    SumWithProduct(this->ReducingFunction, 1., this->ReducingFunction, this->ResolutionTime/3., k);  
    this->EvaluateSourceFunction(k, this->CurrentTime+this->ResolutionTime, this->AuxiliaryFunction); 

    //compute final function
    SumWithProduct(*this->Function, 1., this->ReducingFunction, this->ResolutionTime/6., k);  
    this->CurrentTime += this->ResolutionTime;
}



#endif