#ifndef RUNGEKUTTA_HPP
#define RUNGEKUTTA_HPP

#include <algorithm>
#include <cassert>
#include <functional>
#include <memory>

/*
    Equations implemented here: 
    dy/dt= (A+B)y with y(t0) = y0  

    The time propagation is described by
    exp(f dt)

  f-> source term
*/

template <typename T, typename Scalar_T>
void SumWithProduct(T& Output, const Scalar_T& FirstScalar, const T& FirstAddend, const Scalar_T& SecondScalar, const T& SecondAddend)
{
    assert(Output.end() - Output.begin() == FirstAddend.end() - FirstAddend.begin());
    assert(FirstAddend.end() - FirstAddend.begin() == SecondAddend.end() - SecondAddend.begin());
    
    for(auto OutputIterator=Output.begin(), FirstAddendIterator = FirstAddend.begin(), SecondAddendIterator = SecondAddend.begin();
        (OutputIterator!=Output.end()) && (FirstAddendIterator != FirstAddend.end()) && (SecondAddendIterator != SecondAddend.end());  
          ++OutputIterator, ++FirstAddendIterator, ++SecondAddendIterator ){
                *OutputIterator = FirstScalar*(*FirstAddendIterator) + SecondScalar*(*SecondAddendIterator);
    }
}


//Warning! This class does not own Function, that must be deleted outside!!
template<typename T>
class RungeKutta
{
    private:
        T* Function;            //beginning and end of a step, it has the function at a particular time y(n)
        T AuxiliaryFunction;   //It has the function incremented by k, for the propagator in the 4 steps i.e. yn+h/2*k1
	    T ReducingFunction;    //Adds up all the contributions to the function until the end of the step to get y(n+1) = yn+h/6*(k1+2*k2+2*k3+k4)
        T k;
        double InitialTime = 0.;
        double ResolutionTime = 0.;
        double CurrentTime = 0;
        double FinalTime;
        //LambdaForSourceTerm EvaluateSourceFunction;
        //LambdaForInitialCondition EvaluateInitialCondition;
        std::function<void(T&)> EvaluateInitialCondition;
        std::function<void(T&, const double&, const T&)> EvaluateSourceFunction;
	public:
        RungeKutta(){};

        RungeKutta(const RungeKutta& RK) = default;
        RungeKutta& operator=(const RungeKutta& RK) = default;

        RungeKutta(RungeKutta&& RK) = default;
        RungeKutta& operator=(RungeKutta&& RK) = default;

        RungeKutta(T& Function_, const std::function<void(T&)>& EvaluateInitialCondition_, const std::function<void(T&, const double&, const T&)>& EvaluateSourceFunction_); 
        void initialize(T& Function_, const std::function<void(T&)>& EvaluateInitialCondition_, const std::function<void(T&, const double&, const T&)>& EvaluateSourceFunction_); 
        void Propagate(); 
        const T& get_Function() const; 
        T& get_Function(); 
        const double& get_CurrentTime() const {return CurrentTime;};
        const double& get_ResolutionTime() const {return ResolutionTime;};
        void set_InitialTime(const double& InitialTime_){InitialTime = InitialTime_;}
        void set_ResolutionTime(const double& ResolutionTime_){ResolutionTime = ResolutionTime_;}
};



//template<typename T, typename LambdaForSourceTerm, typename LambdaForInitialCondition>
//RungeKutta<T, LambdaForSourceTerm, LambdaForInitialCondition>::RungeKutta
template<typename T>
RungeKutta<T>::RungeKutta(T& Function_, const std::function<void(T&)>& EvaluateInitialCondition_, 
                        const std::function<void(T&, const double&, const T&)>& EvaluateSourceFunction_)
{
    this->initialize(Function_, EvaluateInitialCondition_, EvaluateSourceFunction_);
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
    //k=f(tn,yn)
    EvaluateSourceFunction(k, CurrentTime, *Function);


    CurrentTime += ResolutionTime;
}


template<typename T>
const T& RungeKutta<T>::get_Function() const 
{
    return *(this->Function);
}

template<typename T>
T& RungeKutta<T>::get_Function() 
{
    return *(this->Function);
}


//template<typename T>
//auto make_RungeKutta(const std::function<void(T&)>& InitialCondition_, const std::function<void(T&, const double&, const T&)>& SourceTerm_) 
//-> RungeKutta<T>
//{
//    return RungeKutta<T>(InitialCondition_, SourceTerm_, 0., 0.);
//}



#endif