
#include "mdContainers/mdContainers.hpp"
#include "MeshGrid/MeshGrid.hpp"



class kGradient {
    private:
        mdarray<double, 1> Weight;    
        mdarray<std::vector<int>, 2> ikpb;
        std::shared_ptr<MeshGrid<k>> kmesh;
        std::vector<std::vector<int>> ikshell;
        int nshells = 0;
    public:
        kGradient(const MeshGrid<k> kmesh__);
        void initialize();

};

int alpha( const size_t& j );
int beta( const size_t& j );
void Calculate_nshellsAndweights(int& nshells, mdarray<double,1>& Weight, 
                                 const MeshGrid<k>& kmesh, const std::vector<std::vector<int>>& ikshell);
auto GradientMatrix(const size_t& nshells, const MeshGrid<k>& kmesh, const std::vector<std::vector<int>>& ik_sorted);
mdarray<std::vector<int>, 2> Find_kpb(const MeshGrid<k>& kmesh, const std::vector<std::vector<int>>& ikshell);