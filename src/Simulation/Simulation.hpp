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

void print_bandstructure(const std::vector<std::vector<double>>& bare_kpath, Operator<std::complex<double>> Hamiltonian);

/// @brief This class contains all the variables that are used in the simulations
class Simulation
{
    private:
    public:
        /// Object containing the main operators for the simulation (H0 and r)
        Material material_;
        /// Object containing all the parameters that needs to be set; they get read or 
        /// simply they get some default values
        std::shared_ptr<Simulation_parameters> ctx_;
        /// Full Hamiltonian: @f$ H_ = H0 + H_{\text{eff}} +E(t)\cdot \Xi @f$, while the term
        /// with the gradient is treated separately
        Operator<std::complex<double>> H_;
        /// This is defined as @f$ V_{nm}(\textbf{k}) = \langle n\textbf{k}| v |m\textbf{k} \rangle @f$
        /// or, in real space, @f$ V_{nm(\textbf{R})} = \langle n\textbf{0} | v | m\textbf{R}\rangle @f$
        std::array<Operator<std::complex<double>>, 3> Velocity_;
        /// Eigenvalues of the Hamiltonian, grouped like: Band_energies[ik](ibnd)
        std::vector<mdarray<double,1>> Band_energies_;
        /// Total laser acting on the system, as a sum over the single lasers
        SetOfLaser setoflaser_;
        /// What we are propagating in time
        Operator<std::complex<double>> DensityMatrix_;
        /// Driver for the time propagation, it defines how we solve the differential equations
        /// and it is fed with our equation
        DESolver<Operator<std::complex<double>>> DEsolver_DM_;
        /// Object to deal with the gradient term. More details in the class
        kGradient kgradient_;
        /// Object to calculate the Coulomb effective Hamiltonian. More details in the class
        Coulomb coulomb_;
        /// Space where we calculate the commutator @f$ [H, \rho] @f$
        Space SpaceOfPropagation_ = k;
        /// Space where we calculate the gradient in k 
        Space SpaceOfPropagation_Gradient_ = R;

        /// Output text file to print the time values where we get the other text files printed
        std::ofstream os_Time_;
        /// Output text file to print the laser in time, electric field      
        std::ofstream os_Laser_;
        /// Output text file to print the laser in time, vector potential
        std::ofstream os_VectorPot_;
        /// Output text file to print the population in time (in bloch gauge) as a sum over all the k points
        std::ofstream os_Pop_;
        /// Output text file to print the mean value of the velocity operator over the state where our system is at time t
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
        double jacobian(const Matrix<double>& A__) const;

        std::string wavelength_or_frequency(const int&);

        template <typename Scalar_T>
        friend void SumWithProduct(Operator<std::complex<double>>& Output__, 
                    const Scalar_T& FirstScalar__, 
                    const Operator<std::complex<double>>& FirstAddend__, 
                    const Scalar_T& SecondScalar__, 
                    const Operator<std::complex<double>>& SecondAddend__);

};

// == template<typename T>
// == std::complex<double> Trace(BlockMatrix<T>& O1__, BlockMatrix<T>& O2__)
// == {
// ==     //calculate trace of product of two operators in R space as \sum_(n,R) (O1(R)O2(-R))_(nn)
// ==     auto& ci = MeshGrid::ConvolutionIndex[{Operator<std::complex<double>>::MeshGrid_Null->get_id(), 
// ==                                               O1__.get_MeshGrid()->get_id(), 
// ==                                               O2__.get_MeshGrid()->get_id()}];
// ==     if(ci.get_Size(0) == 0 ){
// ==         MeshGrid::Calculate_ConvolutionIndex(*(Operator<std::complex<double>>::MeshGrid_Null), 
// ==                                                  *(O1__.get_MeshGrid()), 
// ==                                                  *(O2__.get_MeshGrid()));
// ==     }
// == 
// ==     std::complex<double> Trace = 0.;
// ==     for(int iblock=0; iblock<O1__.get_nblocks(); ++iblock){
// ==         for(int irow=0; irow<O1__.get_nrows(); irow++){
// ==             for(int icol=0; icol<O1__.get_ncols(); icol++){   
// ==                 Trace += O1__[ci(0,iblock)](irow, icol)*O2__[iblock](icol, irow);
// ==             }
// ==         }
// ==     }
// ==     return Trace;
// == }

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
