#include "/home/gcistaro/2negf/header/Geometry/Matrix.hpp"

int main()
{
    //std::cout << "MATRIX A!!\n";
    for(int i=0; i<100; i++){
        int nrows = 300, ncols = 300;
        Matrix<double> A({nrows,ncols});

        std::cout << " SIZE OF A: " << nrows*ncols*sizeof(decltype(A(0,0)));
        std::cout << "MATRIX A DONE!\n";
        Matrix<double> B(ncols,nrows);
        std::cout << "MATRIX B DONE!\n";
        Matrix<double> C;//(ncols, nrows);
        //Matrix_gemm(C,A,B);
	C=A*B;
	A.initialize(3,3);
	A(0,0)=1;
	A(0,1)=1;
	A(0,2)=1;
	A(1,0)=-1;
	A(1,1)=1;
	A(1,2)=-1;
	A(2,0)=1;
	A(2,1)=-1;
	A(2,2)=-1;
	auto invA = A.inverse();
	std::cout << "A: \n";
	std::cout << A; 
	std::cout << "invA: \n";
	std::cout << invA;
	std::cout << "\n" << A*invA;	

    }
}

