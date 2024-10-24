#ifndef DESOLVER_HPP
#define DESOLVER_HPP

#include <algorithm>
#include <cassert>
#include <functional>
#include <memory>

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


															// -- GENERAL FUNCTIONS -- //


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

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


													// -- DIFFERENTIAL EQUATION SOLVER CLASS -- //

template<typename T>
class DESolver{
	//private:
	protected:
		T *Function;                  // Function is declared as a pointer
		double InitialTime = 0.;
        double ResolutionTime = 0.;
        double CurrentTime = 0;
        double FinalTime;
        //LambdaForSourceTerm EvaluateSourceFunction;
        //LambdaForInitialCondition EvaluateInitialCondition;
        std::function<void(T&)> EvaluateInitialCondition;
        std::function<void(T&, const double&, const T&)> EvaluateSourceFunction;

	public:
		DESolver(){};

        DESolver(const DESolver& DES) = default;
        DESolver& operator=(const DESolver& DES) = default;

        DESolver(DESolver&& DES) = default;
        DESolver& operator=(DESolver&& DES) = default;

        // this has to be declared separately in each derived class since it calls initialize, a function defined differently for each method
        //DESolver(T& Function_, const std::function<void(T&)>& EvaluateInitialCondition_, const std::function<void(T&, const double&, const T&)>& EvaluateSourceFunction_); 

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
const T& DESolver<T>::get_Function() const 
{
    return *Function;                         
}

template<typename T>
T& DESolver<T>::get_Function() 
{
    return *Function;
}


//template<typename T>
//auto make_RungeKutta(const std::function<void(T&)>& InitialCondition_, const std::function<void(T&, const double&, const T&)>& SourceTerm_) 
//-> RungeKutta<T>
//{
//    return RungeKutta<T>(InitialCondition_, SourceTerm_, 0., 0.);
//}


#endif