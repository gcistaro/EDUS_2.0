#include <cassert>
#include "ConvertUnits.hpp"
#include "Wannier.hpp"

#include "Operator.hpp"

class Material
{
    private:


    public:
        enum MODEL{wannier} model;

        Operator<std::complex<double>> H;
        std::array< Operator<std::complex<double>>, 3> r;


        Material(){};
        Material(const std::string& Filename)
        {
            model = wannier;
   
            Wannier wannier(Filename);
            //we move the resources from wannier.r and wannier.H, anyway we throw wannier right after.            
            Basis LatticeVectors(Matrix<double>(wannier.UnitCell).transpose());
            Coordinate<R>::add_Basis(LatticeVectors, "LatticeVectors");
            
            Basis ReciprocalLatticeVectors(2.*pi*Matrix<double>(wannier.UnitCell).inverse());
            Coordinate<k>::add_Basis(ReciprocalLatticeVectors, "LatticeVectors");
            
            r[0].get_Operator_R().initialize(wannier.r[0]);
            r[1].get_Operator_R().initialize(wannier.r[1]);
            r[2].get_Operator_R().initialize(wannier.r[2]);
            H.get_Operator_R().initialize(wannier.H);

            Convert_iterable(r[0].get_Operator_R(), Angstrom, Bohr);
            Convert_iterable(r[1].get_Operator_R(), Angstrom, Bohr);
            Convert_iterable(r[2].get_Operator_R(), Angstrom, Bohr);
            Convert_iterable(H.get_Operator_R(), ElectronVolt, Rydberg);

            MeshGrid<R> aux_mg(wannier.Rmesh, "LatticeVectors");

            for(auto& v : aux_mg.get_mesh()){
                std::cout << v << std::endl;
            }
            H.get_Operator_R().set_MeshGrid(aux_mg);
            r[0].get_Operator_R().set_MeshGrid(aux_mg);
            r[1].get_Operator_R().set_MeshGrid(aux_mg);
            r[2].get_Operator_R().set_MeshGrid(aux_mg);
        }

	void print_info()
	{
        std::cout << Coordinate<R>::get_Basis("LatticeVectors");
	    std::cout << Coordinate<k>::get_Basis("LatticeVectors");
	}




};
