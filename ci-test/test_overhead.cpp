#include "initialize.hpp"
#include "MPIindex/MPIindex.hpp"
#include "core/profiler.hpp"
#include "core/print_timing.hpp"
#include "Operator/Operator.hpp"
int main()
{
    initialize();
    MPIindex<3> mpindex;
    BlockMatrix H;
    int i1=100, i2=10, i3=10;
    H.initialize(k, i1,i2,i3);
{
    PROFILE("FIRST LOOP");
    for(int counter = 0; counter < 100; counter++){
        #pragma omp parallel for
        for(int i1_=0; i1_<i1; i1_++) {
            for(int i2_=0; i2_<i2; i2_++) {
                for(int i3_=0; i3_<i3; i3_++) {
                    auto alpha = H(i1_, i2_, i3_);
                }
            }
        }
    }
}
std::cout << "First done.\n";
{
{
    PROFILE("SECOND LOOP");
    for(int counter = 0; counter < 100; counter++){
        #pragma omp parallel for
        for(int i1_=0; i1_<i1; i1_++) {
            for(int i2_=0; i2_<i2; i2_++) {
                for(int i3_=0; i3_<i3; i3_++) {
                    auto alpha = H(i3_+i3*(i2_+i2*i1_));
                }
            }
        }
    }
}

    
}   
    print_timing(1);

}