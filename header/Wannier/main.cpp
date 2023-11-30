#include "Wannier.hpp"
int main()
{
    std::string name = "TBgraphene";

    Wannier wannier(name);
    name = "WhatIRead_tb.dat";
    wannier.Print(name);
}