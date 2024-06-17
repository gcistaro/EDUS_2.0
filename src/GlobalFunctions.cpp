#include "Constants.hpp"
std::string LatticeVectors(const Space& space) 
{ 
    return ( space == k ? "k_LatticeVectors" : "R_LatticeVectors" ); 
};
