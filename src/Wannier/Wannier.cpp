#include "Wannier/Wannier.hpp"
Wannier::Wannier(const std::string& FileName)
{
    PROFILE("Wannier");
    std::stringstream TotalFileName;
    TotalFileName << FileName << "_tb.dat";    
    ParseWannier(TotalFileName.str(), NumberOfBands, NumberOfRpoints,
                  UnitCell, Degeneracy, Rmesh, H, r);
}

void Wannier::Print(const std::string& FileName)
{
    PrintWannier(FileName, NumberOfBands, NumberOfRpoints, 
                UnitCell, Degeneracy, Rmesh, H, r);
}
