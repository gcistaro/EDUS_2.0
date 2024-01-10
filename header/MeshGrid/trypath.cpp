#include "MeshGrid.hpp"

int main()
{
    Matrix<double> M(3,3);
    M(0,0) = 1.;   M(0,1) = 0;   M(0,2) = 0;
    M(1,0) = 0;   M(1,1) = 1.;   M(1,2) = 0;
    M(2,0) = 0;   M(2,1) = 0;   M(2,2) = 1.;
    std::cout <<"M:\n"<< M << std::endl;

    Basis Cartesian;
    Cartesian.initialize(M);
    std::cout << " CAARTESIAN:: \n" << Cartesian;
    Coordinate<k>::add_Basis(Cartesian, "Cartesian");
    M(0,0) = +1.;   M(0,1) = +0.5;   M(0,2) = -1.;
    M(1,0) = -1.;   M(1,1) = +1.5;   M(1,2) = -1.;
    M(2,0) = +1.;   M(2,1) = +0.5;   M(2,2) = +1.;
    Basis LatticeVectors;
    LatticeVectors.initialize(M);
    Coordinate<k>::add_Basis(LatticeVectors, "LatticeVectors");

    std::vector<Coordinate<k>> Path;
    Path.push_back(Coordinate<k>(0,0,0));
    Path.push_back(Coordinate<k>(0,1,0));
    Path.push_back(Coordinate<k>(0,1,2));

    MeshGrid<k> m1(Path, 0.1);
    auto& vectors = m1.get_mesh();
    std::cout << "PATH:::\n";
    for(int i=0; i<vectors.size(); i++){
        std::cout << vectors[i].get("Cartesian");
    }
}