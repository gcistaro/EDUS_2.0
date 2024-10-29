#include <iomanip>
#include "initialize.hpp"
#include "Model/Model.hpp"
#include "core/projectdir.hpp"

int main()
{
    initialize();

    std::stringstream inputfile;
    inputfile << ProjectDirectory <<  "/tb_models/TBgraphene";
    Material model(inputfile.str());

    std::vector<Coordinate> path;

    path.push_back( Coordinate( 0.000, 0.000, 0.000, LatticeVectors(k) ));
    path.push_back( Coordinate( 0.500, 0.500, 0.000, LatticeVectors(k) ));
    path.push_back( Coordinate( 1./3., 2./3., 0.000, LatticeVectors(k) ));
    path.push_back( Coordinate( 0.000, 0.000, 0.000, LatticeVectors(k) ));
    

    MeshGrid MeshGridPath(k, path, 0.1);
    auto& vectors = MeshGridPath.get_mesh();

    for(auto& v : vectors){
        std::cout << v.get(LatticeVectors(k));
    }
    
    model.H.dft(MeshGridPath.get_mesh(),+1,false);
    std::cout << "dft executed.\n";

    std::vector<mdarray<double,1>> Eigenvalues;
    BlockMatrix<std::complex<double>> Eigenvectors;
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


