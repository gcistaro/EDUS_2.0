#include "mdContainers/mdContainers.hpp"
#include "Operator/Operator.hpp"
#include <type_traits>
//#include "Constants.hpp"

template <typename T>
class Decay
{
    private:
        /// Matrix of decays
        T Gamma_;
        /// boolean that is true if Gamma was set
        bool Gamma_set_ = false;
        /// Density matrix for the ground state
        Operator<std::complex<double>> DM0_;
        /// boolean that is true if DM0_ was set
        bool DM0_set_ = false;
        /// Density matrix variation rescaled by decays
        Operator<std::complex<double>> Gamma_DM_DM0_;

    public:

        void set_Gamma(const T& Gamma__);
        void set_DM0(const Operator<std::complex<double>>& DM0__);
        void set_Gamma_DM_DM0(const Operator<std::complex<double>>& DM__);

        const T& get_Gamma() const;
        const Operator<std::complex<double>>& get_DM0() const;
        const Operator<std::complex<double>>& get_Gamma_DM_DM0() const;
};


template <typename T>
void Decay<T>::set_Gamma(const T& Gamma__){
    if (Gamma_set_){
        throw std::runtime_error("Gamma_ has already been set and is immutable.");
    }
    Gamma_ = Gamma__;
    Gamma_set_ = true;
}

template <typename T>
void Decay<T>::set_DM0(const Operator<std::complex<double>>& DM0__){
    if (DM0_set_){
        throw std::runtime_error("DM0_ has already been set and is immutable.");
    }
    DM0_ = DM0__;
    DM0_set_ = true;
}

template <typename T>
void Decay<T>::set_Gamma_DM_DM0(const Operator<std::complex<double>>& DM__){
    auto DMR = DM__.get_Operator(R);
    auto& Gamma_DM_DM0_R = Gamma_DM_DM0_.get_Operator(R);
    Gamma_DM_DM0_R.fill(0.0);
    auto DM0R = DM0_.get_Operator(R);

    for (int iR = 0; iR < DMR.get_nblocks(); iR++){
        for (int irow = 0; irow < DMR.get_nrows(); irow++){
            for (int icol = 0; icol < DMR.get_ncols(); icol++){
                if constexpr (std::is_same<T, mdarray<double, 2>>::value){
                    assert(Gamma_.get_nrows() == DMR.get_nrows() && DMR.get_ncols() == DM0R.get_ncols());
                    assert(Gamma_.get_nrows() == DMR.get_nrows() && DMR.get_ncols() == DM0R.get_ncols());
                    Gamma_DM_DM0_R(iR, irow, icol) = Gamma_(irow, icol) * (DMR(iR, irow, icol) - DM0R(iR, irow, icol));
                }
                else{
                    Gamma_DM_DM0_R(iR, irow, icol) = Gamma_ * (DMR(iR, irow, icol) - DM0R(iR, irow, icol));
                }
            }
        }
    }

    Gamma_DM_DM0_.go_to_wannier();

}

template <typename T>
const T& Decay<T>::get_Gamma() const{
    if (!Gamma_set_){
        throw std::runtime_error("Gamma_ has not been set yet.");
    }
    return Gamma_;
}

template <typename T>
const Operator<std::complex<double>>& Decay<T>::get_DM0() const{
    if (!DM0_set_){
        throw std::runtime_error("DM0_ has not been set yet.");
    }
    return DM0_;
}

template <typename T>
const Operator<std::complex<double>>& Decay<T>::get_Gamma_DM_DM0() const{
    return Gamma_DM_DM0_;
}
