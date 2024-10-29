#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include "ConvertUnits.hpp"
#include "Model/Model.hpp"
#include "RungeKutta/RungeKutta.hpp"
#include "Laser/Laser.hpp"
#include "ostream.hpp"
#include "kGradient/kGradient.hpp"
#include "Coulomb/Coulomb.hpp"
#include "initialize.hpp"
#include "Json/json.hpp"

class Simulation
{
    private:
    public:
        Material material;

        Operator<std::complex<double>> H;
        std::array<Operator<std::complex<double>>, 3> Velocity;
        std::vector<mdarray<double,1>> Band_energies;
        double FermiEnergy = 0.;
        SetOfLaser setoflaser;
        Operator<std::complex<double>> DensityMatrix;
        RungeKutta<Operator<std::complex<double>>> RK_object;
        int PrintResolution; //steps needed to print a variable
        double FinalTime;
        double InitialTime;

        kGradient kgradient;
        Coulomb coulomb;
        Space SpaceOfPropagation = k;
        std::string JsonFile;
        std::string tb_model;
        std::ofstream os_Laser;
        std::ofstream os_VectorPot;
        std::ofstream os_Pop;
        std::ofstream os_Time;
        std::ofstream os_Velocity;
        
        template<class T>
        Simulation(const std::string& FileName, const T& arg_meshinit);
        Simulation(const std::string& JsonFileName);
        void SettingUp_EigenSystem();
        void print_grids();
        void Calculate_TDHamiltonian(const double& time, const bool& erase_H);
        void Calculate_Velocity();
        void Propagate();
        void do_onestep();
        void print_recap();
        bool PrintObservables(const double& time) const;
        void Print_Population();
        void Print_Velocity();


        template <typename Scalar_T>
        friend void SumWithProduct(Operator<std::complex<double>>& Output, 
                    const Scalar_T& FirstScalar, 
                    const Operator<std::complex<double>>& FirstAddend, 
                    const Scalar_T& SecondScalar, 
                    const Operator<std::complex<double>>& SecondAddend);

};

template<typename T>
std::complex<double> Trace(BlockMatrix<T>& O1, BlockMatrix<T>& O2)
{
    //calculate trace of product of two operators in R space as \sum_(n,R) (O1(R)O2(-R))_(nn)
    auto& ci = MeshGrid::ConvolutionIndex[{Operator<std::complex<double>>::MeshGrid_Null->get_id(), 
                                              O1.get_MeshGrid()->get_id(), 
                                              O2.get_MeshGrid()->get_id()}];
    if(ci.get_Size(0) == 0 ){
        MeshGrid::Calculate_ConvolutionIndex(*(Operator<std::complex<double>>::MeshGrid_Null), 
                                                 *(O1.get_MeshGrid()), 
                                                 *(O2.get_MeshGrid()));
    }

    std::complex<double> Trace = 0.;
    for(int iblock=0; iblock<O1.get_nblocks(); ++iblock){
        for(int irow=0; irow<O1.get_nrows(); irow++){
            for(int icol=0; icol<O1.get_ncols(); icol++){   
                Trace += O1[ci(0,iblock)](irow, icol)*O2[iblock](icol, irow);
            }
        }
    }
    return Trace;
}


template<typename T>
std::vector<std::complex<double>> TraceK(BlockMatrix<T>& O__)
{
    //calculate trace over k P(n) = \sum_k O(k)_{nn}
    std::vector<std::complex<double>> TraceK_(O__.get_nrows(), 0.);
    for(int iblock=0; iblock<O__.get_nblocks(); ++iblock){
        for(int ibnd=0; ibnd<O__.get_nrows(); ++ibnd){
            TraceK_[ibnd] += O__[iblock](ibnd, ibnd);
        }
    }
#ifdef EDUS_MPI
    std::vector<std::complex<double>> TraceK_reduced(O__.get_nrows(), 0.);
    kpool_comm.reduce(&TraceK_[0], &TraceK_reduced[0], TraceK_.size(), MPI_SUM, 0);
    return TraceK_reduced;
#else
    return TraceK_;
#endif
}


std::string wavelength_or_frequency(const nlohmann::json& data);


#include "Simulation_definitions.hpp"

#endif