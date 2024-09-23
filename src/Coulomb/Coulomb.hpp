
#include "mdContainers/mdContainers.hpp"
#include "MeshGrid/MeshGrid.hpp"
#include "Operator/Operator.hpp"


class Coulomb 
{
    private:
        mdarray<std::complex<double>, 6> W; //contains Coulomb Interaction in Wannier basis
        std::shared_ptr<MeshGrid> Rgrid;
        std::shared_ptr<Operator<std::complex<double>>> DM0;
        mdarray<std::complex<double>, 3> HF;
    public:
        Coulomb(){};
        
        Coulomb(const int& nbnd, const std::shared_ptr<MeshGrid>& Rgrid__);
        void initialize(const int& nbnd, const std::shared_ptr<MeshGrid>& Rgrid__);
        Coulomb(const int& nbnd, const std::shared_ptr<MeshGrid>& Rgrid__, const std::array<Operator<std::complex<double>>,3>& r);
        void initialize(const int& nbnd, const std::shared_ptr<MeshGrid>& Rgrid__, const std::array<Operator<std::complex<double>>,3>& r);

        void EffectiveHamiltonian(BlockMatrix<std::complex<double>>& H_, const BlockMatrix<std::complex<double>>& DM,
                                  const bool& EraseH_ );        
};



