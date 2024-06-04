#include "Geometry/Matrix.hpp"
#include "Geometry/Vector.hpp"
#include "core/print_timing.hpp"

#define watch(x)  (#x)
constexpr std::complex<double> im(0.,1.);
int main()
{
    PROFILE_START("NEGF");
    std::cout << "Testing matrix multiplication...\n";
    double threshold = 1.e-14;
    Matrix<double> A(3,3);
    Matrix<double> B(3,3);
    Matrix<double> C(3,3);
    Matrix<double> D(3,4);
    Matrix<double> dummy_C(3,3);
    Matrix<double> dummy_E(3,4);

    
    C.fill(0);
    dummy_C.fill(0);

    A(0,0) = 2.;     A(0,1) = 3.;     A(0,2) = 1.;
    A(1,0) = 3.;     A(1,1) = 3.;     A(1,2) = 1.;
    A(2,0) = 2.;     A(2,1) = 4.;     A(2,2) = 1.;

    B(0,0) = 1.;     B(0,1) = 0.;     B(0,2) = 1.;
    B(1,0) = 0.;     B(1,1) = 2.;     B(1,2) = 1.;
    B(2,0) = 1.;     B(2,1) = 1.;     B(2,2) = 1.;

    D(0,0) = 1.;     D(0,1) = 0.;     D(0,2) = 1.;      D(0,2) = 4.;
    D(1,0) = 0.;     D(1,1) = 1.;     D(1,2) = 0.;      D(1,2) = 2.;
    D(2,0) = 1.;     D(2,1) = 1.;     D(2,2) = 1.;      D(2,2) = 1.;
    C = A*B;
    
    auto E = A*D;
    for(int irow=0; irow<dummy_C.get_nrows(); irow++){
        for(int icol=0; icol<dummy_C.get_ncols(); icol++){
            for(int index=0; index<A.get_ncols(); index++){
                dummy_C(irow, icol) += A(irow, index)*B(index, icol);
            }
        }
    }

    std::cout <<"E rows: " << E.get_nrows() << " E cols: " << E.get_ncols() <<"\n";
    for(int irow=0; irow<dummy_E.get_nrows(); irow++){
        for(int icol=0; icol<dummy_E.get_ncols(); icol++){
            for(int index=0; index<A.get_ncols(); index++){
                dummy_E(irow, icol) += A(irow, index)*D(index, icol);
            }
        }
    }

    //--multiplication of matrix of same size
    std::cout << "C:\n" << C << std::endl;
    std::cout << "dummy_C:\n" << dummy_C << std::endl;
    std::cout << (C-dummy_C).norm() << std::endl;

    if((C-dummy_C).norm() > threshold){
        exit(1);
    }

    //--multiplication of matrix of different size
    std::cout << "E:\n" << E << std::endl;
    std::cout << "dummy_E:\n" << dummy_E << std::endl;
    std::cout << (E-dummy_E).norm() << std::endl;

    if((E-dummy_E).norm() > threshold){
        exit(1);
    }



    std::cout << "Testing inverse...\n";
    C = A.inverse();

    dummy_C(0,0) = -1.;     dummy_C(0,1) = +1.;     dummy_C(0,2) = +0.;
    dummy_C(1,0) = -1.;     dummy_C(1,1) = +0.;     dummy_C(1,2) = +1.;
    dummy_C(2,0) = +6.;     dummy_C(2,1) = -2.;     dummy_C(2,2) = -3.;

    if((C-dummy_C).norm() > threshold){
        exit(1);
    }

    B(0,0) = 1.;     B(0,1) = 0.;     B(0,2) = 0.;
    B(1,0) = 0.;     B(1,1) = 1.;     B(1,2) = 0.;
    B(2,0) = 0.;     B(2,1) = 0.;     B(2,2) = 1.;
    
    if((A*C-B).norm() >threshold || (C*A-B).norm() > threshold){
        exit(1);
    }


    //test Matrix-vector multiplication
    Vector<double> v; 
    v.initialize_n(1,2,3);

    std::cout << "B\n" << B;
    std::cout << "v\n" << v;
    std::cout <<  "B*v\n" << B*v;

    PROFILE_STOP("NEGF");



    //--test matrix matrx multiplication with dimension 1 
    A = Matrix<double>(1,6);
    B = Matrix<double>(6,1);
    A(0,0) = 2.;
    A(0,1) = 2.;
    A(0,2) = 0.;
    A(0,3) = 0.; 
    A(0,4) = 0.;
    A(0,5) = 0.;

    B(0,0) = .25; B(1,0) = .25; B(2,0) = 0.; B(3,0) = 0.; B(4,0) = 0.; B(5,0) = 0.;
    std::cout <<"A*B: "<< A*B << std::endl;
    print_timing(2);
}
