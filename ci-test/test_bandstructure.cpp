#include <iomanip>
#include "Model/Model.hpp"
int main()
{

    Material model("/home/gcistaro/NEGF/tb_models/TBgraphene");

    std::vector<Coordinate<k>> path;

    path.push_back(Coordinate<k>(0.000, 0.000, 0.000, "LatticeVectors"));
    path.push_back(Coordinate<k>(0.500, 0.500, 0.000, "LatticeVectors"));
    path.push_back(Coordinate<k>(1./3., 2./3., 0.000, "LatticeVectors"));
    path.push_back(Coordinate<k>(0.000, 0.000, 0.000, "LatticeVectors"));

    MeshGrid<k> MeshGridPath(path, 0.1);
    auto& vectors = MeshGridPath.get_mesh();

    for(auto& v : vectors){
        std::cout << v.get("LatticeVectors");
    }
    
    std::cout << model.H.get_Operator_R()[0];
    model.H.dft(MeshGridPath.get_mesh(),+1);

    std::vector<mdarray<double,1>> Eigenvalues;
    BlockMatrix<std::complex<double>, k> Eigenvectors;
    model.H.get_Operator_k().diagonalize(Eigenvalues, Eigenvectors);

    std::ofstream Output;
    Output.open("BANDSTRUCTURE.txt");
    Output << "#k-number    energy(eV)\n";
    for(int ik=0; ik<Eigenvalues.size(); ik++){
        for(int iband=0; iband<Eigenvalues[ik].get_Size(0); ++iband){
            Output << std::setw(6) << ik;
            Output << std::setw(15) << std::setprecision(6) << Convert(Eigenvalues[ik](iband),AuEnergy,ElectronVolt) << std::endl;
        }
    }
    Output.close();

}


