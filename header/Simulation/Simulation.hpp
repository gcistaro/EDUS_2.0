#include "Model/Model.hpp"
#include "RungeKutta/RungeKutta.hpp"
#include "Laser/Laser.hpp"

class Simulation
{
    private:
        Material material;

        Operator<std::complex<double>> H;
        std::vector<mdarray<double,1>> Band_energies;
        double FermiEnergy;
        Laser laser;
        Operator<std::complex<double>> DensityMatrix;
        RungeKutta<BlockMatrix<std::complex<double>, R>> RK_object;

    public:
        Simulation(const std::string& FileName, const double& Radius)
        {
            material = Material(FileName);
            std::cout << material.r[0].get_Operator_R()[0] << std::endl;
            laser.set_Amplitude(10);
            laser.set_InitialTime(0.);
            laser.set_Omega(0.1);
            laser.set_NumberOfCycles(3);
            laser.set_Polarization(Coordinate<k>(1,0,0));

            //auto& DensityMatrix = RK_object.get_Function();
            //setting up grids
            auto MasterRgrid = std::make_shared<MeshGrid<R>>(MeshGrid<R>(Radius));//spherical grid for exponentially decaying DM

            //get U and Udagger from Hamiltonian
            DensityMatrix.initialize_fft(*MasterRgrid, material.H.get_Operator_R().get_nrows());
            Operator<std::complex<double>>::temp_k = DensityMatrix.get_Operator_k();
            H = material.H;//to initialize dimensions.
            SettingUp_EigenSystem();
            
            auto& Uk = Operator<std::complex<double>>::EigenVectors;
            auto& UkDagger = Operator<std::complex<double>>::EigenVectors_dagger;
            std::function<void(BlockMatrix<std::complex<double>, R>&)> InitialCondition = [&](BlockMatrix<std::complex<double>, R>& DM)
            {
                //DM.Operator_k.initialize(Uk.get_nblocks(), Uk.get_nrows(), Uk.get_ncols());
                DM = DensityMatrix.get_Operator_R();
                DensityMatrix.Operator_k.fill(0.);
                //filling matrix in Bloch gauge
                for(int ik=0; ik<Uk.get_nblocks(); ++ik){
                    for(int iband=0; iband<Uk.get_nrows(); iband++){
                        if(this->Band_energies[ik](iband) < FermiEnergy-threshold){
                            DensityMatrix.Operator_k(ik, iband, iband) = 1.;
                        }
                    }
                }
                DensityMatrix.lock_gauge(bloch);
                DensityMatrix.lock_space(k);
                DensityMatrix.go_to_wannier();
                DensityMatrix.go_to_R();
            };

            std::function<void(BlockMatrix<std::complex<double>, R>&, double const&, BlockMatrix<std::complex<double>, R> const&)> SourceTerm = 
            [&](BlockMatrix<std::complex<double>, R>& Output, const double& time, const BlockMatrix<std::complex<double>, R>& Input)
            {
                Output.fill(0.*im);
                Calculate_TDHamiltonian(time);
                convolution(Output, -im, H.get_Operator_R(), Input);
                convolution(Output, +im, Input, H.get_Operator_R());
            };

            //RK_object.set_Function(DensityMatrix.get_Operator_R());
            RK_object.initialize(InitialCondition, SourceTerm);

            RK_object.set_InitialTime(0.);
            RK_object.set_ResolutionTime(0.01);

            RK_object.Propagate();
        };

        void SettingUp_EigenSystem()
        {
            //diagonalizing Hk on the k vectors of the grid of the Density matrix.
            auto& MasterkGrid = DensityMatrix.get_Operator_k().get_MeshGrid()->get_mesh();
            material.H.dft(MasterkGrid, +1);
            
            material.H.get_Operator_k().diagonalize(Band_energies, Operator<std::complex<double>>::EigenVectors);

            auto& Uk = Operator<std::complex<double>>::EigenVectors;
            auto& UkDagger = Operator<std::complex<double>>::EigenVectors_dagger;
            //HermitianTranspose(Operator<std::complex<double>>::EigenVectors, 
            //                   Operator<std::complex<double>>::EigenVectors_dagger);

            //doing conjugate transpose using mkl...
            UkDagger.initialize(Uk.get_nblocks(), Uk.get_ncols(), Uk.get_nrows());
            //no stride
            //size_t lda = Uk.get_nrows();
            //size_t ldb = Uk.get_ncols();
            //mkl_zomatcopy('R', 'C', Uk.get_nrows(), Uk.get_ncols(), std::complex<double>(1.), 
            //              &Uk(0,0,0), lda, &UkDagger(0,0,0), ldb); 
            for(int ik=0; ik<UkDagger.get_nblocks(); ++ik){
                for(int ir=0; ir<UkDagger.get_nrows(); ++ir){
                    for(int ic=0; ic<UkDagger.get_ncols(); ++ic){
                        UkDagger[ik](ir, ic) = conj(Uk[ik](ic,ir));
                    }
                }
            }

            for(int ik=0; ik<UkDagger.get_nblocks(); ++ik){
                std::cout << UkDagger[ik]*Uk[ik];
            }

        };



        void print_grids()
        {
            auto& MasterRgrid = DensityMatrix.get_Operator_R().get_MeshGrid();
            auto& kgrid = DensityMatrix.get_Operator_k().get_MeshGrid();
            auto Rfft = std::make_shared<MeshGrid<R>>(fftPair<k, R>(*kgrid));//Rspace as fourier space of kgrid
            std::ofstream os;
            os.open("MasterR.txt");
            for (auto &R : MasterRgrid->get_mesh())
            {
                os << R.get("Cartesian");
            }
            os.close();
            os.open("k.txt");
            for (auto &k : kgrid->get_mesh())
            {
                os << k.get("Cartesian");
            }
            os.close();
            auto MeshRfft = fftPair<k, R>(*kgrid);
            os.open("R.txt");
            for (auto &R : Rfft->get_mesh())
            {
                os << R.get("Cartesian");
            }
            os.close();
        }
        

        void Calculate_TDHamiltonian(const double& time)
        {
            auto& HR = H.get_Operator_R();
            auto& H0R = material.H.get_Operator_R();
            auto& xR = material.r[0].get_Operator_R();
            auto& yR = material.r[1].get_Operator_R();
            auto& zR = material.r[2].get_Operator_R();
            std::cout << xR[0];
            auto& las = laser(time).get("Cartesian");
            for(int iblock=0; iblock<H.get_Operator_R().get_nblocks(); ++iblock){
                for(int irow=0; irow<H.get_Operator_R().get_nrows(); ++irow){
                    for(int icol=0; icol<H.get_Operator_R().get_ncols(); ++icol){
                        HR(iblock, irow, icol) = H0R(iblock, irow, icol);
                        HR(iblock, irow, icol) += las[0]*material.r[0].get_Operator_R()(iblock, irow, icol);
                        HR(iblock, irow, icol) += las[1]*material.r[1].get_Operator_R()(iblock, irow, icol);
                        HR(iblock, irow, icol) += las[2]*material.r[2].get_Operator_R()(iblock, irow, icol);
                    }
                }
            }
        };
};

