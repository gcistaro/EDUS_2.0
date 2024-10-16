// ask mikhail about initialization

#ifndef RUNGEKUTTA_HPP
#define RUNGEKUTTA_HPP

#include <algorithm>
#include <cassert>
#include <functional>
#include <memory>

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


//Warning! This class does not own Function, that must be deleted outside!!
template<typename T>
class RungeKutta
{
    private:
        T *Function;       			
        std::array<T,4> fns;
        int ifn1 = 0, ifn2 = 1, ifn3 = 2, ifn4 = 3;
        std::array<double,4> beta = {55/24, -59/24, 37/24, -3/8};
	    
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
        void set_EvaluateInitialCondition(const std::function<void(T&)>& EvaluateInitialCondition_);

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

    fns[ifn1] = *Function, fns[ifn2] = *Function, fns[ifn3] = *Function, fns[ifn4] = *Function;

    std::fill(fns[ifn1].begin(), fns[ifn1].end(), 0.);
    std::fill(fns[ifn2].begin(), fns[ifn2].end(), 0.);
    std::fill(fns[ifn3].begin(), fns[ifn3].end(), 0.);
    std::fill(fns[ifn4].begin(), fns[ifn4].end(), 0.);

}

template<typename T>
void RungeKutta<T>::set_EvaluateInitialCondition(const std::function<void(T&)>& EvaluateInitialCondition_)
{
    EvaluateInitialCondition = EvaluateInitialCondition_;
}


template<typename T>
void RungeKutta<T>::Propagate()
{

    // y_n = y_{n-1} + h*b_1f(t_{n-1}, y_{n-1})
    SumWithProduct(*Function, 1., *Function, ResolutionTime*beta[ifn1], fns[ifn1]);

    // y_{n} += h*b_2f(t_{n-2}, y_{n-2})
    SumWithProduct(*Function, 1., *Function, ResolutionTime*beta[ifn2], fns[ifn2]);

    // y_{n} += h*b_3f(t_{n-3}, y_{n-3})
    SumWithProduct(*Function, 1., *Function, ResolutionTime*beta[ifn3], fns[ifn3]);

    // y_{n} += h*b_4f(t_{n-4}, y_{n-4})
    SumWithProduct(*Function, 1., *Function, ResolutionTime*beta[ifn4], fns[ifn4]);

    ifn4 = ifn3;
    ifn3 = ifn2;
    ifn2 = ifn1;
    ifn1 = 6 - (ifn2 + ifn3 + ifn4);			// ifn1 + ifn2 + ifn3 + ifn4 = 6

    EvaluateSourceFunction(fns[ifn1], CurrentTime, *Function);

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