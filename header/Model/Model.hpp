#include <cassert>
#include "ConvertUnits.hpp"
#include "Wannier.hpp"
#include "MeshGrid.hpp"
#include "BlockMatrix.hpp"

class Model
{
    private:


    public:
        enum MODEL{wannier} model;
        std::array<BlockMatrix<std::complex<double>,R>, 3> r;
        BlockMatrix<std::complex<double>,R> H;
        MeshGrid<R> Master_MeshGrid;


        Model(){};
        Model(const std::string& Filename)
        {
            model = wannier;
   
            Wannier wannier(Filename);
            //we move the resources from wannier.r and wannier.H, anyway we throw wannier right after.            
            Basis LatticeVectors(Matrix<double>(wannier.UnitCell).transpose());
            Coordinate<R>::add_Basis(LatticeVectors, "LatticeVectors");
            r[0].initialize(wannier.r[0]);
            r[1].initialize(wannier.r[1]);
            r[2].initialize(wannier.r[2]);
            H.initialize(wannier.H);

            Convert_iterable(r[0], Angstrom, Bohr);
            Convert_iterable(r[1], Angstrom, Bohr);
            Convert_iterable(r[2], Angstrom, Bohr);
            Convert_iterable(H, ElectronVolt, Rydberg);

            MeshGrid<R> aux_mg(wannier.Rmesh, "LatticeVectors");
            
            H.set_MeshGrid(aux_mg);
            r[0].set_MeshGrid(aux_mg);
            r[1].set_MeshGrid(aux_mg);
            r[2].set_MeshGrid(aux_mg);

            //Master_MeshGrid.initialize(10.);//TODO: change 10. with a variable we can read
            //std::cout << "Number of grids: " << MeshGrid<R>::get_counter_id() << std::endl;
        }


};
