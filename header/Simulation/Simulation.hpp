#include "Model/Model.hpp"
#include "RungeKutta/RungeKutta.hpp"


template<typename LambdaForSourceTerm, typename LambdaForInitialCondition>
class Simulation
{
    private:
        Material material;
        //Operator<std::complex<double>> DensityMatrix;
        RungeKutta<Operator<std::complex<double>>, LambdaForSourceTerm, LambdaForInitialCondition> RK_object;

    public:
        Simulation(const LambdaForInitialCondition& EvaluateInitialCondition_, 
                   const LambdaForSourceTerm& EvaluateSourceFunction_) 
        : RK_object(make_RungeKutta<Operator<std::complex<double>>>(EvaluateInitialCondition_, EvaluateSourceFunction_))
        {
            material = Material("/home/gcistaro/2negf_github/tb_models/TBgraphene");

            auto& DensityMatrix = RK_object.get_Function();
            //setting up grids
            auto& MasterRgrid = DensityMatrix.get_Operator_R().get_MeshGrid();
            MasterRgrid = std::make_shared<MeshGrid<R>>(MeshGrid<R>(3.));//spherical grid for exponentially decaying DM
            
            auto& kgrid = DensityMatrix.get_Operator_k().get_MeshGrid();
            kgrid = std::make_shared<MeshGrid<k>>(std::move(fftPair<R, k>(*MasterRgrid)));//fftPair for k space

            auto Rfft = std::make_shared<MeshGrid<R>>(fftPair<k, R>(*kgrid));//Rspace as fourier space of kgrid
            print_grids();   
            material.print_info();
        };





        void print_grids()
        {
            auto& DensityMatrix = RK_object.get_Function();
            auto& MasterRgrid = DensityMatrix.get_Operator_R().get_MeshGrid();
            auto& kgrid = DensityMatrix.get_Operator_k().get_MeshGrid();
            auto Rfft = std::make_shared<MeshGrid<R>>(fftPair<k, R>(*kgrid));//Rspace as fourier space of kgrid
            std::ofstream os;
            os.open("MasterR.txt");
            for (auto &R : MasterRgrid->get_mesh())
            {
                os << R.get("Cartesian");
            }
            os.close();
            os.open("k.txt");
            for (auto &k : kgrid->get_mesh())
            {
                os << k.get("Cartesian");
            }
            os.close();
            auto MeshRfft = fftPair<k, R>(*kgrid);
            os.open("R.txt");
            for (auto &R : Rfft->get_mesh())
            {
                os << R.get("Cartesian");
            }
            os.close();
        }
        
};



template<typename T, typename LambdaForSourceTerm, typename LambdaForInitialCondition>
auto make_Simulation(const LambdaForInitialCondition& InitialCondition_, const LambdaForSourceTerm& SourceTerm_) 
-> Simulation<LambdaForSourceTerm, LambdaForInitialCondition>
{
    return Simulation<LambdaForSourceTerm, LambdaForInitialCondition>(InitialCondition_, SourceTerm_);
}
