#ifndef DESOLVER_HPP
#define DESOLVER_HPP

#include <algorithm>
#include <cassert>
#include <functional>
#include <memory>
#include "Operator/Operator.hpp"

template <typename T, typename Scalar_T>
void SumWithProduct(T& Output, const Scalar_T& FirstScalar, const T& FirstAddend, const Scalar_T& SecondScalar, const T& SecondAddend)
{
    assert(Output.end() - Output.begin() == FirstAddend.end() - FirstAddend.begin());
    assert(FirstAddend.end() - FirstAddend.begin() == SecondAddend.end() - SecondAddend.begin());
  
    #pragma omp parallel for
    for( int i=0; i<Output.end()-Output.begin(); ++i) {
        *(Output.begin()+i) = FirstScalar*(*(FirstAddend.begin()+i)) + SecondScalar*(*(SecondAddend.begin()+i));
    }
}


template<typename T>
class DESolver{
    protected:
		T* Function;                  //It contains the function at a particular time y(n)
		double InitialTime = 0.;
        double ResolutionTime = 0.;
        double CurrentTime = 0;
        double FinalTime;
        std::function<void(T&)> EvaluateInitialCondition;
        std::function<void(T&, const double&, const T&)> EvaluateSourceFunction;

        std::array<T,4> aux_Function;                        //temporary arrays to store quantities needed in DE numerical methods 

        std::array<int, 4> index = {0, 1, 2, 3};              // to rotate indices in AB method, if not we need to do many copies
        std::array<double,4> beta = {55./24., -59./24., 37./24., -3./8.}; //coefficients in AB method

        void Propagate_RK4();
        void Propagate_AB4();
        SolverType type;

    public:
        DESolver(){};

        DESolver(const DESolver& DEsolver__) = default;
        DESolver& operator=(const DESolver& DEsolver__) = default;

        DESolver(DESolver&& DEsolver__) = default;
        DESolver& operator=(DESolver&& DEsolver__) = default; 


        DESolver(T& Function_, const std::function<void(T&)>& EvaluateInitialCondition__, const std::function<void(T&, const double&, const T&)>& EvaluateSourceFunction__, SolverType type__); 

        const T& get_Function() const; 
        T& get_Function(); 
        const double& get_CurrentTime() const {return CurrentTime;};
        const double& get_ResolutionTime() const {return ResolutionTime;};
        void set_InitialTime(const double& InitialTime_){InitialTime = InitialTime_;}
        void set_ResolutionTime(const double& ResolutionTime_){ResolutionTime = ResolutionTime_;}
        void set_EvaluateInitialCondition(const std::function<void(T&)>& EvaluateInitialCondition_);
        void initialize(T& Function_, const std::function<void(T&)>& EvaluateInitialCondition_, const std::function<void(T&, const double&, const T&)>& EvaluateSourceFunction_, SolverType type__);
        void Propagate();
        void set_aux_Function(const T& a, const T& b, const T& c, const T& d){aux_Function[0] = a, aux_Function[1] = b, aux_Function[2] = c, aux_Function[3] = d;}
        void set_type(SolverType t){type = t;}
        SolverType get_type(){return type;}	
};

template<typename T>
void DESolver<T>::initialize(T& Function__, const std::function<void(T&)>& EvaluateInitialCondition__, const std::function<void(T&, const double&, const T&)>& EvaluateSourceFunction__, SolverType type__)
{
    std::cout << "Not specialized\n"; exit(0);
    type = type__;
    Function = &Function__;
    EvaluateInitialCondition = EvaluateInitialCondition__;
    EvaluateSourceFunction = EvaluateSourceFunction__;
    EvaluateInitialCondition(*Function);

    for( auto& aux_ : aux_Function) {
        aux_ = *Function;
        std::fill(aux_.begin(), aux_.end(), 0.);
    }
}

template<typename T>
DESolver<T>::DESolver(T& Function__, const std::function<void(T&)>& EvaluateInitialCondition__, 
                        const std::function<void(T&, const double&, const T&)>& EvaluateSourceFunction__, SolverType type__)
{
    initialize(Function__, EvaluateInitialCondition__, EvaluateSourceFunction__, type__);
}


template<typename T>
void DESolver<T>::Propagate_RK4()
{
    /*
        Equations implemented here: 
        dy/dt= f(t,y) with y(t0) = y0  
      In RK4:
      y(n+1) = y(n)+h/6*(k1+2*k2+2*k3+k4)
      k1=f(tn,yn)
      k2=f(tn+h/2, yn+h/2*k1)
      k3=f(tn+h/2, yn+h/2*k2)
      k4=f(tn+h, yn+h*k3)
      the calculation of k- is defined as EvaluateSourceTerm.

      f-> source term    */

    auto& k = aux_Function[0];                  //k1, k2, k3 in RK method
    auto& AuxiliaryFunction = aux_Function[1];  //temporary values of function, second argument of EvaluateSourceFunction
    auto& ReducingFunction = aux_Function[2];   //sum up contributions from a step to next

    //k1=f(tn,yn)
    EvaluateSourceFunction(k, CurrentTime, *Function);

    //k2=f(tn+h/2,yn+h/2*k1)
    SumWithProduct(AuxiliaryFunction, 1., *Function, ResolutionTime/2., k);     
    SumWithProduct(ReducingFunction, 1., *Function, ResolutionTime/6., k);     
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

template<typename T>
void DESolver<T>::Propagate_AB4()
{
    /*
        Equations implemented here: 
        dy/dt= f(t,y) with y(t0) = y0  
      In AB4: 

      y(n)  = y(n-1) + h* 55./24.*f(t(n-1), y(n-1)) +
                     - h* 59./24.*f(t(n-2), y(n-2)) +
                     + h* 37./24.*f(t(n-3), y(n-3)) +
                     - h*  9./24.*f(t(n-4), y(n-4))

      f-> source term    
      aux_Function is used to store the previous steps propagators    
      index rotates in order to avoid copies                          */

    //get f(t(n-1), y(n-1))
    EvaluateSourceFunction(aux_Function[index[0]], CurrentTime, *Function);

    // y_n = y_{n-1} + h*b_1*f(t_{n-1}, y_{n-1})
    SumWithProduct(*Function, 1., *Function, ResolutionTime*beta[0], aux_Function[index[0]]);

    // y_{n} += h*b_2*f(t_{n-2}, y_{n-2})
    SumWithProduct(*Function, 1., *Function, ResolutionTime*beta[1], aux_Function[index[1]]);

    // y_{n} += h*b_3*f(t_{n-3}, y_{n-3})
    SumWithProduct(*Function, 1., *Function, ResolutionTime*beta[2], aux_Function[index[2]]);

    // y_{n} += h*b_4*f(t_{n-4}, y_{n-4})
    SumWithProduct(*Function, 1., *Function, ResolutionTime*beta[3], aux_Function[index[3]]);

    index[3] = index[2];
    index[2] = index[1];
    index[1] = index[0];
    index[0] = 6 - (index[1] + index[2] + index[3]);


    CurrentTime += ResolutionTime;
}

template<typename T>
void DESolver<T>::Propagate(){
    if (type == RK4){Propagate_RK4();}
    else if (type == AB4){Propagate_AB4();}
}

template<typename T>
const T& DESolver<T>::get_Function() const 
{
    return *Function;                         
}

template<typename T>
T& DESolver<T>::get_Function() 
{
    return *Function;
}



template<>
void DESolver<Operator<std::complex<double>>>::initialize(Operator<std::complex<double>>& Function_, 
                const std::function<void(Operator<std::complex<double>>&)>& EvaluateInitialCondition_, 
                const std::function<void(Operator<std::complex<double>>&, 
                const double&, const Operator<std::complex<double>>&)>& EvaluateSourceFunction_,
                SolverType type__);


#endif