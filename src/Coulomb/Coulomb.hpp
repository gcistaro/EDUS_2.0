
#include "mdContainers/mdContainers.hpp"
#include "MeshGrid/MeshGrid.hpp"
#include "RytovaKeldysh/RytovaKeldysh.hpp"

class Coulomb 
{
    private:
        mdarray<std::complex<double>, 6> W; //contains Coulomb Interaction in Wannier basis
        std::shared_ptr<MeshGrid> Rgrid;
        RytovaKeldysh RytKel;
    public:
        Coulomb(){};
        
        Coulomb(const int& nbnd, const MeshGrid& Rgrid__);
        void initialize(const int& nbnd, const MeshGrid& Rgrid__);

        void EffectiveHamiltonian(BlockMatrix<std::complex<double>>& H_, const BlockMatrix<std::complex<double>>& DM,
                                  const bool& EraseH_ );        
};



