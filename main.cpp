//to install:::::::
//g++ main.cpp -L/opt/intel/oneapi/mkl/2022.2.1/lib/intel64 -lmkl_rt


#include <iostream> 
#include <vector>
#include <complex>
#include <memory>   //for smart pointers
#include <cassert>  //for assertions, avoiding errors for inconsistencies
#include <iomanip>  //for setprecision and setw
#include <time.h> //for clock

#include "mkl.h"//for blas-lapack

double c1,c2,c3;
bool direct, cartesian;
//#include "header/Matrix.hpp"
//#include "header/Geometry.hpp"
//#include "header/MeshGrid.hpp"
//#include "header/Tensor.hpp"


#define watch(x)  (#x)
constexpr std::complex<double> im(0.,1.);
int main()
{
    int NumberOfK = 1;
    int NumberOfBands = 10;
    l2::Tensor3 Hamiltonian(NumberOfK, NumberOfBands);
    //Operator Position(NumberOfK, NumberOfBands);
    //Operator DensityMatrix(NumberOfK, NumberOfBands);
    l2::Tensor3 Eigenvectors(NumberOfK, NumberOfBands);

    Hamiltonian[0](0,0) = 1.; Hamiltonian[0](0,1) = 3.; 
    Hamiltonian[0](1,0) = 2.; Hamiltonian[0](1,1) = 1.; 
    Matrix<std::complex<double>> B(NumberOfBands);
    B(0,0) = 2.; B(0,1) = 1.; 
    B(1,0) = .5; B(1,1) = 5.; 
 



    clock_t t =clock();
    Matrix<std::complex<double>> C;  
    Matrix<std::complex<double>> A(NumberOfBands);  
    C=A*B;
    C=Hamiltonian[0]*B;
    t = clock() - t;    
    printf ("* took me %f seconds).\n",((float)t)/CLOCKS_PER_SEC);
    //std::cout << C;


    //std::cout << "test of r3" << std::endl;
    //r3::Basis A;
    //Matrix<double> A_(3);
    //A_(0,0)= 1.05;     A_(0,1)= 2.07;     A_(0,2)= 0.0;  
    //A_(1,0)= 0.99;     A_(1,1)=-2.15;     A_(1,2)= 1.27;
    //A_(2,0)= 0.00;     A_(2,1)= 0.00;     A_(2,2)= 7.5;
    //A.initialize(A_);
    //std::cout << A;
    ////test of r3 is Basis is working. Checked with Mathematica.
    //r3::Vector::DirectLatticeVectors.initialize(A_);
//
    //r3::Vector v(c1=0,c2=1,c3=0,direct=true,cartesian=true);
    //auto c=v.get_cart();
    //std::cout << c[0] << " " << c[1] << " " << c[2] << std::endl;
    //c=v.get_notcart();
    //std::cout << c[0] << " " << c[1] << " " << c[2] << std::endl;
//
//
    //std::array<int,3> S={10,10,1};
    //MeshGrid R(S,direct=true);
    //
    //std::vector<r3::Vector> Rbas; 
    //Rbas.push_back(r3::Vector(c1= 1,c2= 0,c3=0,direct=true,cartesian=false));
    //Rbas.push_back(r3::Vector(c1=-1,c2= 0,c3=0,direct=true,cartesian=false));
    //Rbas.push_back(r3::Vector(c1= 0,c2= 1,c3=0,direct=true,cartesian=false));
    //Rbas.push_back(r3::Vector(c1= 0,c2=-1,c3=0,direct=true,cartesian=false));
    //for(auto& R1 : R.get_mesh()){
    //    for(auto& R2 : Rbas){
    //        auto R_ = R1+R2;
    //        int index = R.where(R_);
    //    }
    //}
}
