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
    Coordinate<R>::add_Basis(Cartesian, "Cartesian");
    M(0,0) = +1.;   M(0,1) = +0.5;   M(0,2) = -1.;
    M(1,0) = -1.;   M(1,1) = +1.5;   M(1,2) = -1.;
    M(2,0) = +1.;   M(2,1) = +0.5;   M(2,2) = +1.;
    Basis LatticeVectors;
    LatticeVectors.initialize(M);
    Coordinate<R>::add_Basis(LatticeVectors, "LatticeVectors");
   
    //sphere mesh testing:
    MeshGrid<R> m1(double(40.));
 /* 
    for(int iR1=0; iR1<m1.get_TotalSize(); iR1++){
        std::cout << "m1[iR1]:\n" << m1[iR1].get("LatticeVectors");
        std::cout << "m1[iR1].norm(): " << m1[iR1].norm() << std::endl;
    }
*/
    MeshGrid<R>::Calculate_ConvolutionIndex(m1,m1,m1);
    auto ci = MeshGrid<R>::get_ConvolutionIndex(m1,m1,m1);

   //std::cout << "id: " << m1.get_id() << std::endl;
/*
    //square mesh testing:

    MeshGrid<R> m1(std::array<int,3>{3,2,3});
    MeshGrid<R>::Calculate_ConvolutionIndex(m1,m1,m1);
    auto ci = MeshGrid<R>::get_ConvolutionIndex(m1,m1,m1);
*/    
/*    
    std::cout << std::endl;
    for(int iR1=0; iR1<m1.get_TotalSize(); iR1++){
        for(int iR3=0; iR3<m1.get_TotalSize(); iR3++){
                std::cout << "iR1 << << iR3: " <<iR1 << " " << " " << iR3 << std::endl;

                std::cout << m1[iR1].get("LatticeVectors") << m1[iR3].get("LatticeVectors");
                std::cout << "m1[iR1]-m3[iR3].norm(): ";
                std::cout << std::setprecision(15) << (m1[iR1]-m1[iR3]).norm() << std::endl;
                std::cout << "ci(iR1,iR3): " << ci(iR1,iR3) << std::endl;
            std::cout << std::endl;
            std::cout << std::endl;
            std::cout << std::endl;
            std::cout << std::endl;
       }
    }
    */
}
