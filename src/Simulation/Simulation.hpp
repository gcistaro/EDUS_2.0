#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include "ConvertUnits.hpp"
#include "Model/Model.hpp"
#include "DESolver/DESolver.hpp"
#include "Laser/Laser.hpp"
#include "ostream.hpp"
#include "kGradient/kGradient.hpp"
#include "Coulomb/Coulomb.hpp"
#include "initialize.hpp"
#include "Json/json.hpp"
#include "InputVariables/simulation_parameters.hpp"

class Simulation
{
    private:
    public:
        Material material_;
        
        std::shared_ptr<Simulation_parameters> ctx_;
        
        Operator<std::complex<double>> H_;
        std::array<Operator<std::complex<double>>, 3> Velocity_;
        std::vector<mdarray<double,1>> Band_energies_;
        SetOfLaser setoflaser_;
        Operator<std::complex<double>> DensityMatrix_;
        DESolver<Operator<std::complex<double>>> DEsolver_DM_;

        kGradient kgradient_;
        Coulomb coulomb_;
        Space SpaceOfPropagation_ = k;
        std::ofstream os_Laser_;
        std::ofstream os_VectorPot_;
        std::ofstream os_Pop_;
        std::ofstream os_Time_;
        std::ofstream os_Velocity_;
        
        Simulation(){};
        Simulation(std::shared_ptr<Simulation_parameters>& ctx__);
        void SettingUp_EigenSystem();
        void print_grids();
        void Calculate_TDHamiltonian(const double& time__, const bool& erase_H__);
        void Calculate_Velocity();
        void Propagate();
        void do_onestep();
        void print_recap();
        bool PrintObservables(const double& time__, const bool& use_sparse = true);
        void Print_Population();
        void Print_Velocity();
        int get_it(const double& time__) const;
        int get_it_sparse(const double& time__) const;

        std::string wavelength_or_frequency(const int&);
        
        template <typename Scalar_T>
        friend void SumWithProduct(Operator<std::complex<double>>& Output__, 
                    const Scalar_T& FirstScalar__, 
                    const Operator<std::complex<double>>& FirstAddend__, 
                    const Scalar_T& SecondScalar__, 
                    const Operator<std::complex<double>>& SecondAddend__);

};

template<typename T>
std::complex<double> Trace(BlockMatrix<T>& O1__, BlockMatrix<T>& O2__)
{
    //calculate trace of product of two operators in R space as \sum_(n,R) (O1(R)O2(-R))_(nn)
    auto& ci = MeshGrid::ConvolutionIndex[{Operator<std::complex<double>>::MeshGrid_Null->get_id(), 
                                              O1__.get_MeshGrid()->get_id(), 
                                              O2__.get_MeshGrid()->get_id()}];
    if(ci.get_Size(0) == 0 ){
        MeshGrid::Calculate_ConvolutionIndex(*(Operator<std::complex<double>>::MeshGrid_Null), 
                                                 *(O1__.get_MeshGrid()), 
                                                 *(O2__.get_MeshGrid()));
    }

    std::complex<double> Trace = 0.;
    for(int iblock=0; iblock<O1__.get_nblocks(); ++iblock){
        for(int irow=0; irow<O1__.get_nrows(); irow++){
            for(int icol=0; icol<O1__.get_ncols(); icol++){   
                Trace += O1__[ci(0,iblock)](irow, icol)*O2__[iblock](icol, irow);
            }
        }
    }
    return Trace;
}


template<typename T>
std::vector<std::complex<double>> TraceK(BlockMatrix<T>& O__)
{
    //calculate trace over k P(n) = \sum_k O(k)_{nn}
    std::vector<std::complex<double>> TraceK(O__.get_nrows(), 0.);
    for(int iblock=0; iblock<O__.get_nblocks(); ++iblock){
        for(int ibnd=0; ibnd<O__.get_nrows(); ++ibnd){
            TraceK[ibnd] += O__[iblock](ibnd, ibnd);
        }
    }
#ifdef EDUS_MPI
    std::vector<std::complex<double>> TraceK_reduced(O__.get_nrows(), 0.);
    kpool_comm.reduce(&TraceK[0], &TraceK_reduced[0], TraceK.size(), MPI_SUM, 0);
    return TraceK_reduced;
#else
    return TraceK;
#endif
}




#include "Simulation_definitions.hpp"

#endif