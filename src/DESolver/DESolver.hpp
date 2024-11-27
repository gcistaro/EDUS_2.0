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

        // RK4 functions
        T AuxiliaryFunction;   
        T ReducingFunction;    
        T k;

        // AB4 functions
        std::array<T,4> fns;
        int ifn1 = 0, ifn2 = 1, ifn3 = 2, ifn4 = 3;
        std::array<double,4> beta = {55./24., -59./24., 37./24., -3./8.};
        
	public:
        enum Type {RK4, AB4};
        Type type;
        
		DESolver(){};

        DESolver(const DESolver& DES) = default;
        DESolver& operator=(const DESolver& DES) = default;

        DESolver(DESolver&& DES) = default;
        DESolver& operator=(DESolver&& DES) = default; 


        DESolver(T& Function_, const std::function<void(T&)>& EvaluateInitialCondition_, const std::function<void(T&, const double&, const T&)>& EvaluateSourceFunction_, Type type); 

        // this has to be declared separately in each derived class since it calls initialize, a function defined differently for each method
        //DESolver(T& Function_, const std::function<void(T&)>& EvaluateInitialCondition_, const std::function<void(T&, const double&, const T&)>& EvaluateSourceFunction_); 

        const T& get_Function() const; 
        T& get_Function(); 
        const double& get_CurrentTime() const {return CurrentTime;};
        const double& get_ResolutionTime() const {return ResolutionTime;};
        void set_InitialTime(const double& InitialTime_){InitialTime = InitialTime_;}
        void set_ResolutionTime(const double& ResolutionTime_){ResolutionTime = ResolutionTime_;}
        void set_EvaluateInitialCondition(const std::function<void(T&)>& EvaluateInitialCondition_);
        void initialize(T& Function_, const std::function<void(T&)>& EvaluateInitialCondition_, const std::function<void(T&, const double&, const T&)>& EvaluateSourceFunction_); 
        void Propagate();
        void Propagate_RK4();
        void Propagate_AB4();
        void setFnsAB4(const T& a, const T& b, const T& c, const T& d){fns[0] = a, fns[1] = b, fns[2] = c, fns[3] = d;}
        void setType(Type t){type = t;}
        Type getType(){return type;}

	
};

template<typename T>
void DESolver<T>::initialize(T& Function_, const std::function<void(T&)>& EvaluateInitialCondition_, const std::function<void(T&, const double&, const T&)>& EvaluateSourceFunction_)
{
    if (type == RK4){
        Function = &Function_;
        EvaluateInitialCondition = EvaluateInitialCondition_;
        EvaluateSourceFunction = EvaluateSourceFunction_;
        EvaluateInitialCondition(*Function);
        AuxiliaryFunction = *Function;
        ReducingFunction = *Function;
        k = *Function; 
    }
    else if (type == AB4){
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
}

template<typename T>
DESolver<T>::DESolver(T& Function_, const std::function<void(T&)>& EvaluateInitialCondition_, 
                        const std::function<void(T&, const double&, const T&)>& EvaluateSourceFunction_, Type type)
{
    setType(type);
    initialize(Function_, EvaluateInitialCondition_, EvaluateSourceFunction_);
}

//template<typename T, typename LambdaForSourceTerm, typename LambdaForInitialCondition>
//RungeKutta<T, LambdaForSourceTerm, LambdaForInitialCondition>::RungeKutta

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


													// -- RungeKutta4 Propagation -- //

template<typename T>
void DESolver<T>::Propagate_RK4(){

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

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


													// -- AdamsBashforth4 Propagation -- //

template<typename T>
void DESolver<T>::Propagate_AB4(){

    // uncomment this for the density matrix. it is commented for testing

    /*
    fns[ifn1] = *Function; fns[ifn2] = *Function; fns[ifn3] = *Function; fns[ifn4] = *Function;
    std::fill(fns[ifn1].begin(), fns[ifn1].end(), 0.);
    std::fill(fns[ifn2].begin(), fns[ifn2].end(), 0.);
    std::fill(fns[ifn3].begin(), fns[ifn3].end(), 0.);
    std::fill(fns[ifn4].begin(), fns[ifn4].end(), 0.);
    */


    /*

    we want to switch the indexes so that in the next step it saves as f(t_{n-1}, y_{n-1}) what is f(t_n, y_n) the step before

    this means that, for instance f(t_{n-4}, y_{n-4}) is what was only 3 steps ago in the previous step, thus

        . f(t_{n-4}, y_{n-4}) -> f(t_{n-3}, y_{n-3})
        . f(t_{n-3}, y_{n-3}) -> f(t_{n-2}, y_{n-2})
        . f(t_{n-3}, y_{n-3}) -> f(t_{n-1}, y_{n-1})

    let us do 2 steps to see that this makes sense

        y_{n+1} = y_n + h*(b1*f(t_{n-1}, y_{n-1}) + b2*f(t_{n-2}, y_{n-2}) + b3*f(t_{n-3}, y_{n-3}) + b4*f(t_{n-4}, y_{n-4}))

                = y_n + h*(b1*fns[0] + b2*fns[1] + b3*fns[2] + b4*fns[3])


        y_{n+2} = y_{n+1} + h*(b1*g(t_{n-1}, y_{n-1}) + b2*g(t_{n-2}, y_{n-2}) + b3*g(t_{n-3}, y_{n-3}) + b4*g(t_{n-4}, y_{n-4}))

                = y_{n+1} + h*(b1*f(t_{n-4}, y_{n-4}) + b2*f(t_{n-1}, y_{n-1}) + b3*f(t_{n-2}, y_{n-2}) + b4*f(t_{n-3}, y_{n-3}))

                = y_{n+1} + h*(b1*fns[3] + b2*fns[0] + b3*fns[1] + b4*fns[2])

    where g(t_{n-1}, y_{n-1}) = f(t_{n-4}, y_{n-4}) = fns[3] needs to be reevaluated in this last step.

    this tells us that

            y_{n+1}    ->   y_{n+2}    ->   ...


            ifn1 = 0        ifn1 = 3
            ifn2 = 1        ifn2 = 0
            ifn3 = 2        ifn3 = 1
            ifn4 = 3        ifn4 = 2

    or in another way, in the next step the indexes change as

            ifn4 = ifn3
            ifn3 = ifn2
            ifn2 = ifn1
            ifn1 = 6 - (ifn2 + ifn3 + ifn4)

    where 6 is the sum of the indexes. in order

    */

    //std::cout << "ifn1 = " << ifn1 << " ifn2 = " << ifn2 << " ifn3 = " << ifn3 << " ifn4 = " << ifn4 << "\n"; 

    // y_n = y_{n-1} + h*b_1*f(t_{n-1}, y_{n-1})
    SumWithProduct(*Function, 1., *Function, ResolutionTime*beta[0], fns[ifn1]);
    //SumWithProduct(*Function, 1., *Function, ResolutionTime*beta[0], fns[0]);

    // y_{n} += h*b_2*f(t_{n-2}, y_{n-2})
    SumWithProduct(*Function, 1., *Function, ResolutionTime*beta[1], fns[ifn2]);
    //SumWithProduct(*Function, 1., *Function, ResolutionTime*beta[1], fns[1]);

    // y_{n} += h*b_3*f(t_{n-3}, y_{n-3})
    SumWithProduct(*Function, 1., *Function, ResolutionTime*beta[2], fns[ifn3]);
    //SumWithProduct(*Function, 1., *Function, ResolutionTime*beta[2], fns[2]);

    // y_{n} += h*b_4*f(t_{n-4}, y_{n-4})
    SumWithProduct(*Function, 1., *Function, ResolutionTime*beta[3], fns[ifn4]);
    //SumWithProduct(*Function, 1., *Function, ResolutionTime*beta[3], fns[3]);

    ifn4 = ifn3;
    ifn3 = ifn2;
    ifn2 = ifn1;
    ifn1 = 6 - (ifn2 + ifn3 + ifn4);

    EvaluateSourceFunction(fns[ifn1], CurrentTime, *Function);

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




//template<typename T>
//auto make_RungeKutta(const std::function<void(T&)>& InitialCondition_, const std::function<void(T&, const double&, const T&)>& SourceTerm_) 
//-> RungeKutta<T>
//{
//    return RungeKutta<T>(InitialCondition_, SourceTerm_, 0., 0.);
//}


#endif