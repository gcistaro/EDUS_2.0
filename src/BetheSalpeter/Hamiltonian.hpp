#include "Simulation/Simulation.hpp"

#ifndef BSE_HAM_HPP
#define BSE_HAM_HPP

class BS_Hamiltonian 
{
    private:
        Simulation simulation_;
        Matrix<std::complex<double>> H_BSE_;
        std::shared_ptr<Simulation_parameters> ctx_;
        MultiIndex<3> TwoBodyIndex_;
        //DavidsonSolver davidson_;

    public:
        void initialize(std::shared_ptr<Simulation_parameters> ctx__, const int& iq__);
};

#endif