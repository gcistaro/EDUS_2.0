#include "Geometry/Matrix.hpp"
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
    Matrix<double> dummy_C(3,3);

    C.fill(0);
    dummy_C.fill(0);

    A(0,0) = 2.;     A(0,1) = 3.;     A(0,2) = 1.;
    A(1,0) = 3.;     A(1,1) = 3.;     A(1,2) = 1.;
    A(2,0) = 2.;     A(2,1) = 4.;     A(2,2) = 1.;

    B(0,0) = 1.;     B(0,1) = 0.;     B(0,2) = 1.;
    B(1,0) = 0.;     B(1,1) = 2.;     B(1,2) = 1.;
    B(2,0) = 1.;     B(2,1) = 1.;     B(2,2) = 1.;

    C = A*B;

    for(int irow=0; irow<dummy_C.get_nrows(); irow++){
        for(int icol=0; icol<dummy_C.get_ncols(); icol++){
            for(int index=0; index<A.get_ncols(); index++){
                dummy_C(irow, icol) += A(irow, index)*B(index, icol);
            }
        }
    }

    std::cout << "C:\n" << C << std::endl;
    std::cout << "dummy_C:\n" << dummy_C << std::endl;
    std::cout << (C-dummy_C).norm() << std::endl;

    if((C-dummy_C).norm() > threshold){
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

    PROFILE_STOP("NEGF");
    print_timing(2);
}
