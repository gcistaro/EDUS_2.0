#include "Simulation/Simulation.hpp"

class ScreenedPotential : public Simulation
{
    private: 
        Operator<std::complex<double>> W_;              //screened interaction
        Operator<std::complex<double>> V_;              //bare interaction
        Operator<std::complex<double>> X_;              //it's the response function Chi
        Operator<std::complex<double>> Epsilon_;        //permittivity
        Operator<std::complex<double>> InverseEpsilon_; //inverse of permittivity

        std::string nodename_;

    public: 
        using Simulation::Simulation;

        void initialize();

        void ResponseFunction();
        void Epsilon();
        void BareCoulomb();
        void Calculate();
};
