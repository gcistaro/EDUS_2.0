//g++ -fconcepts -std=c++17 -I../header/ test_readWannier.cpp -o test_wannier.x -L/opt/intel/oneapi/mkl/2022.2.1/lib/intel64 -lmkl_rt

#include "ReadWannier.hpp"
#include "PrintWannier.hpp"
int main()
{
    std::string FileName = "TBgraphene_tb.dat";
    int NumberOfBands;
    int NumberOfRpoints;
    std::vector<int> Degeneracy;
    l2::Tensor3 H;
    std::array<l2::Tensor3,3> r;
    std::cout << "Parsing wannier file...\n";
    ParseWannier(FileName, NumberOfBands, NumberOfRpoints,
                      Degeneracy, H, r);
    PrintWannier("Test_code.dat", NumberOfBands, NumberOfRpoints,
                      Degeneracy, H, r);

}
