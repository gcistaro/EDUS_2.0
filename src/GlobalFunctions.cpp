#include "GlobalFunctions.hpp"

std::string LatticeVectors(const Space& space) 
{ 
    return ( space == k ? "k_LatticeVectors" : "R_LatticeVectors" ); 
};


bool
file_exists(std::string file_name)
{
    std::ifstream ifs(file_name.c_str());
    return ifs.is_open();
}
