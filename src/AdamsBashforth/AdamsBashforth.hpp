#ifndef ADAMSBASHFORTH_HPP
#define ADAMSBASHFORTH_HPP
#include "DESolver/DESolver.hpp"

                                                        // -- ADAMS BASHFORTH DERIVED CLASS -- //

template<typename T>
class AdamsBashforth: public DESolver<T>{
    private:                
        std::array<T,4> fns;
        int ifn1 = 0, ifn2 = 1, ifn3 = 2, ifn4 = 3;
        std::array<double,4> beta = {55./24., -59./24., 37./24., -3./8.};

        using DESolver<T>::Function;
        using DESolver<T>::ResolutionTime;
        using DESolver<T>::CurrentTime;
        using DESolver<T>::EvaluateInitialCondition;
        using DESolver<T>::EvaluateSourceFunction;

    public:
        AdamsBashforth(T& Function_, const std::function<void(T&)>& EvaluateInitialCondition_, const std::function<void(T&, const double&, const T&)>& EvaluateSourceFunction_) 
        : DESolver<T>(Function_, EvaluateInitialCondition_, EvaluateSourceFunction_){};

        void initialize(T& Function_, const std::function<void(T&)>& EvaluateInitialCondition_, const std::function<void(T&, const double&, const T&)>& EvaluateSourceFunction_) override; 
        void Propagate(); 
        void set_fns(const T& a, const T& b, const T& c, const T& d);
};

template<typename T>
void AdamsBashforth<T>::initialize(T& Function_, const std::function<void(T&)>& EvaluateInitialCondition_, const std::function<void(T&, const double&, const T&)>& EvaluateSourceFunction_)
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
void AdamsBashforth<T>::set_fns(const T& a, const T& b, const T& c, const T& d){

    fns[0] = a, fns[1] = b, fns[2] = c, fns[3] = d;
}


template<typename T>
void AdamsBashforth<T>::Propagate()
{

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



#endif