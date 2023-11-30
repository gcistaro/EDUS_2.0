#include "Matrix.hpp"

int main()
{
    for(int i=0; i<100; i++){
        int nrows = 5, ncols = 5;
        Matrix<double> A({nrows,ncols});
        Matrix<double> B(ncols,nrows);

        A(0,0) = ;  A(0,1) = ;  A(0,2) = ;  A(0,3) = ;  A(0,4) = ;
        A(1,0) = ;  A(1,1) = ;  A(1,2) = ;  A(1,3) = ;  A(1,4) = ;
        A(2,0) = ;  A(2,1) = ;  A(2,2) = ;  A(2,3) = ;  A(2,4) = ;
        A(3,0) = ;  A(3,1) = ;  A(3,2) = ;  A(3,3) = ;  A(3,4) = ;
        A(4,0) = ;  A(4,1) = ;  A(4,2) = ;  A(4,3) = ;  A(4,4) = ;

        Matrix<double> C;//(ncols, nrows);
    	C=A*B;
    }
}