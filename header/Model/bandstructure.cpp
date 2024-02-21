#include "Model.hpp"
int main()
{

    std::cout << "Reading wannier file...\n";
    Material model("TBgraphene");
    std::cout << "read!\n";
    auto& HR = model.H.get_Operator_R();

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

    model.H.get_Operator_k().test_submatrix();
    std::vector<mdarray<double,1>> Eigenvalues;
    BlockMatrix<std::complex<double>, k> Eigenvectors;
    for(int ik=0; ik<MeshGridPath.get_mesh().size(); ik++){
    std::cout <<" model.H.get_Operator_k():: \n\n\n" <<  model.H.get_Operator_k().Values << std::endl;
    }
    diagonalize(model.H.get_Operator_k(), Eigenvalues, Eigenvectors);

    std::ofstream Output;
    Output.open("BANDSTRUCTURE.txt");
    for(int ik=0; ik<Eigenvalues.size(); ik++){
        for(int iband=0; iband<Eigenvalues[ik].get_Size(0); ++iband){
            Output << ik << " " << Eigenvalues[ik](iband) << std::endl;
        }
    }
    Output.close();

}


