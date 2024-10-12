
#include "mdContainers/mdContainers.hpp"
#include "MeshGrid/MeshGrid.hpp"
#include "Operator/Operator.hpp"
#include "StreamFile.hpp"

class Coulomb 
{
    private:
        mdarray<std::complex<double>, 6> W; //contains Coulomb Interaction in Wannier basis
        std::shared_ptr<MeshGrid> Rgrid;
        Operator<std::complex<double>> DM0;
        mdarray<std::complex<double>, 3> HF;
        bool DoCoulomb = false;
    public:
        Coulomb(){};
        
        Coulomb(const int& nbnd, const std::shared_ptr<MeshGrid>& Rgrid__);
        void initialize(const int& nbnd, const std::shared_ptr<MeshGrid>& Rgrid__);
        Coulomb(const int& nbnd, const std::shared_ptr<MeshGrid>& Rgrid__, const std::array<Operator<std::complex<double>>,3>& r);
        void initialize(const int& nbnd, const std::shared_ptr<MeshGrid>& Rgrid__, const std::array<Operator<std::complex<double>>,3>& r);

        void EffectiveHamiltonian(Operator<std::complex<double>>& H__, const Operator<std::complex<double>>& DM__,
                                  const bool& EraseH__);     

        void set_DM0( const Operator<std::complex<double>>& DM0__ );
        void set_DoCoulomb(const bool& DoCoulomb__);
        const bool& get_DoCoulomb() const;
        bool& get_DoCoulomb();

};



