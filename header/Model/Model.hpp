#include <cassert>
#include "ConvertUnits.hpp"
#include "core/profiler.hpp"
#include "Wannier/Wannier.hpp"

#include "Operator/Operator.hpp"

class Material
{
    private:


    public:
        enum MODEL{wannierTB} model;

        Operator<std::complex<double>> H;
        std::array< Operator<std::complex<double>>, 3> r;


        Material(){};
        Material(const Material& m) = default;
        Material& operator=(const Material& m) = default;
        
        Material(Material&& m) = default;
        Material& operator=(Material&& m) = default;

        Material(const std::string& Filename)
        {
            PROFILE("Material");
            model = wannierTB;
   
            Wannier wannier_(Filename);
            //we move the resources from wannier.r and wannier.H, anyway we throw wannier right after.   
            auto UC = Matrix<double>(wannier_.UnitCell).transpose();
            Convert_iterable(UC, Angstrom, AuLength);   
            Basis LatticeVectors(UC);
            Coordinate<R>::add_Basis(LatticeVectors, "LatticeVectors");
            
            Basis ReciprocalLatticeVectors(2.*pi*Matrix<double>(wannier_.UnitCell).inverse());
            Coordinate<k>::add_Basis(ReciprocalLatticeVectors, "LatticeVectors");
            
            r[0].get_Operator_R().initialize(wannier_.r[0]);
            r[1].get_Operator_R().initialize(wannier_.r[1]);
            r[2].get_Operator_R().initialize(wannier_.r[2]);
            H.get_Operator_R().initialize(wannier_.H);

            Convert_iterable(r[0].get_Operator_R(), Angstrom, AuLength);
            Convert_iterable(r[1].get_Operator_R(), Angstrom, AuLength);
            Convert_iterable(r[2].get_Operator_R(), Angstrom, AuLength);
            Convert_iterable(H.get_Operator_R(), ElectronVolt, AuEnergy);

            MeshGrid<R> aux_mg(wannier_.Rmesh, "LatticeVectors");

            H.get_Operator_R().set_MeshGrid(aux_mg);
            r[0].get_Operator_R().set_MeshGrid(aux_mg);
            r[1].get_Operator_R().set_MeshGrid(aux_mg);
            r[2].get_Operator_R().set_MeshGrid(aux_mg);

            r[0].lock_gauge(wannier);      r[0].lock_space(R);        
            r[1].lock_gauge(wannier);      r[1].lock_space(R);        
            r[2].lock_gauge(wannier);      r[2].lock_space(R);        
            H.lock_gauge(wannier);         H.lock_space(R);         



        }

	void print_info()
	{
        std::cout << Coordinate<R>::get_Basis("LatticeVectors");
	    std::cout << Coordinate<k>::get_Basis("LatticeVectors");
	}




};
