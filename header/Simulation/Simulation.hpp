#include "Model/Model.hpp"
#include "RungeKutta/RungeKutta.hpp"
#include "Laser/Laser.hpp"


class Simulation
{
    private:
        Material material;

        Operator<std::complex<double>> H;
        std::vector<mdarray<double,1>> Band_energies;
        double FermiEnergy;
        Laser laser;
        Operator<std::complex<double>> DensityMatrix;
        //RungeKutta<BlockMatrix<std::complex<double>, k>> RK_object;
        RungeKutta<BlockMatrix<std::complex<double>, R>> RK_object;
        std::ofstream OutLaser;
        std::ofstream OutPos;
    public:
        Simulation(const std::string& FileName, const double& Radius);
        void SettingUp_EigenSystem();
        void print_grids();
        void Calculate_TDHamiltonian(const double& time);
        void Propagate();
        void print_recap();
};

