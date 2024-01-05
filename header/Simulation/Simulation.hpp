#include "Model.hpp"

class Simulation
{
    private:
        Material material;
        Operator<std::complex<double>> DensityMatrix;

    public:
        Simulation()
        {
            material = Material("TBgraphene");

            auto MasterRgrid = DensityMatrix.get_Operator_R().get_MeshGrid();
            MasterRgrid = std::make_shared<MeshGrid<R>>(MeshGrid<R>(3.));

            auto kgrid = DensityMatrix.get_Operator_k().get_MeshGrid();
            kgrid = std::make_shared<MeshGrid<k>>(std::move(fftPair<R, k>(*MasterRgrid)));

            auto Rfft = std::make_shared<MeshGrid<R>>(fftPair<k, R>(*kgrid));

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
            material.print_info();

        }
};
