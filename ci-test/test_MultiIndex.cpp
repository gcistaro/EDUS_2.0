#include <iostream>
#include "MultiIndex/MultiIndex.hpp"

int main()
{
    int a = 1;
    std::array<int,3> Dim = {100,200,300};
    MultiIndex mi(Dim);
    std::cout << mi.oneDindex(98,199,1);
    auto ndindex = mi.nDindex(100);
    std::cout << ndindex[0] << " " << ndindex[1] << " " << ndindex[2] << std::endl;
}