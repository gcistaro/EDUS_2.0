#include "Model/Model.hpp"
#include "RungeKutta/RungeKutta.hpp"
#include "Laser/Laser.hpp"



class Simulation
{
    private:
        Material material;

        Operator<std::complex<double>> H;
        std::array<Operator<std::complex<double>>, 3> Velocity;
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
        void Calculate_Velocity();
        void Propagate();
        void print_recap();
};
template<typename T>
std::complex<double> Trace(BlockMatrix<T, R>& O1, BlockMatrix<T, R>& O2)
{
    auto Minus_MG = std::make_shared<MeshGrid<R>>(Opposite(*(O1.get_MeshGrid())));
    auto& ci = MeshGrid<R>::ConvolutionIndex1[{Operator<std::complex<double>>::MeshGrid_Null->get_id(), 
                                              O2.get_MeshGrid()->get_id(), 
                                              Minus_MG->get_id()}];
    if(ci.get_Size(0) == 0 ){
        MeshGrid<R>::Calculate_ConvolutionIndex1(*(Operator<std::complex<double>>::MeshGrid_Null), *(O2.get_MeshGrid()), *Minus_MG);
    }

    std::complex<double> Trace = 0.;
    for(int iblock=0; iblock<O1.get_nblocks(); ++iblock){
        for(int irow=0; irow<O1.get_nrows(); irow++){
            for(int icol=0; icol<O1.get_ncols(); icol++){   
                Trace += O1[iblock](irow, icol)*O2[ci(0,iblock)](icol, irow);
            }
        }
    }
    return Trace;
}
