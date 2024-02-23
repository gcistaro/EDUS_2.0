

#ifndef MESHGRID_HPP
#define MESHGRID_HPP

#include <vector>
#include <algorithm>

#include "../Constants.hpp"
#include "../Geometry/Coordinate.hpp"


struct ShellInfo
{
    double Radius;
    int NumberOfPoints;
    int StartingIndex;
};


enum TypeMeshGrid{cube, read_, sphere, path};
 
template<Space space>
class MeshGrid{
    private:
        std::vector<Coordinate<space>> mesh;
        
        //cube parameter
        std::array<int,3> Size;
        
        //sphere parameter
        double maxRadius = 0;
        std::vector<ShellInfo> shellinfo;

        //auxiliary vector for out of bounds (empty vector)
        Coordinate<space> EmptyVector;

        size_t TotalSize=0;
        TypeMeshGrid type;
        size_t id; //each MeshGrid has a unique id,  to exchange the indices with others
        static size_t counter_id; //neded to count objects
    public:
        static std::map<std::array<size_t,3>, mdarray<int,2> > ConvolutionIndex;//to call it: [{id1,id2,id3}][{iR1,iR3}]

        MeshGrid(){};
        MeshGrid(const MeshGrid& m);
        MeshGrid<space>& operator=(const MeshGrid<space>& mg);

        MeshGrid(MeshGrid&& m);
        MeshGrid<space>& operator=(MeshGrid<space>&& mg);

        MeshGrid(const std::vector<Coordinate<space>>& ReadMesh); 
        MeshGrid(const mdarray<double,2>& bare_mg, const std::string& KeyForBasis);
        MeshGrid(const std::array<int,3>& Size_);
        MeshGrid(const double& Radius_);
        MeshGrid(const std::vector<Coordinate<space>>& PathPoints, const double& resolution);

        void initialize(const std::array<int,3>& Size_);
        void initialize(const double& Radius_);
        void initialize(const mdarray<double,2>& bare_mg, const std::string& KeyForBasis);
        void initialize(const std::vector<Coordinate<space>>& PathPoints, const double& resolution);

        std::array<double,3> get_absmax() const;
        
        template<Space space1, Space space2>
        friend MeshGrid<space2> fftPair(const MeshGrid<space1>& KnownMG);
        
        Coordinate<space>& operator[](const int& i);
        const Coordinate<space>& operator[](const int& i) const;

        int find(const Coordinate<space>& v) const;
        Coordinate<space> reduce(const Coordinate<space>& v) const;

        const std::vector<Coordinate<space>>& get_mesh() const;
        const std::array<int,3>& get_Size() const;
        size_t get_TotalSize() const;
        size_t get_id() const;
        static size_t get_counter_id(){return counter_id;}; //neded to count objects

        static void Calculate_ConvolutionIndex(const MeshGrid& m1, const MeshGrid& m2, const MeshGrid& m3);
        static mdarray<int,2>& get_ConvolutionIndex(const MeshGrid& m1, const MeshGrid& m2, const MeshGrid& m3){
            return ConvolutionIndex[{m1.get_id(), m2.get_id(), m3.get_id()}];
        }
};

#include "MeshGrid_definitions.hpp"

#endif
