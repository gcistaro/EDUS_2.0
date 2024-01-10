#include "Model.hpp"

int main()
{
    Material model("TBgraphene");
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
    
}


