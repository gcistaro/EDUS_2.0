#ifndef KGRADIENT_HPP
#define KGRADIENT_HPP

#include "mdContainers/mdContainers.hpp"
#include "MeshGrid/MeshGrid.hpp"
#include "MPIindex/MPIindex.hpp"
#include "core/mpi/Communicator.hpp"


class kGradient {
    private:
        mdarray<double, 1> Weight;    
        std::vector<std::vector<std::vector<int>>> ikpb;
        std::shared_ptr<MeshGrid> kmesh;
        std::shared_ptr<MeshGrid> Rmesh;
        std::vector<std::vector<int>> ikshell;
        int nshells = 0;

        std::vector<std::vector<std::vector<int>>> Find_kpb(const MeshGrid& kmesh, const std::vector<std::vector<int>>& ikshell);

        //---mpi
        MPIindex<3> mpindex;

    public:
        kGradient(){};
        kGradient(const MeshGrid& kmesh__);
        
        void initialize(const MeshGrid& kmesh__);
        void initialize();

        template<typename T, typename U>
        void Calculate(T& DerivativeFunction, const T& Function, 
                          const U& direction, const bool& EraseOutput) const;
};

int alpha( const size_t& j );
int beta( const size_t& j );
void Calculate_nshellsAndweights(int& nshells, mdarray<double,1>& Weight, 
                                 const MeshGrid& kmesh, const std::vector<std::vector<int>>& ikshell);
Matrix<double> GradientMatrix(const size_t& nshells, const MeshGrid& kmesh, const std::vector<std::vector<int>>& ik_sorted);
std::vector<std::vector<int>> SortInShells(const MeshGrid& kmesh);

#include "kGradient_definitions.hpp"

#endif