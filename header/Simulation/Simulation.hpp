#include "Model/Model.hpp"
#include "RungeKutta/RungeKutta.hpp"


class Simulation
{
    private:
        Material material;

        std::vector<mdarray<double,1>> Band_energies;
        double FermiEnergy;
        //Operator<std::complex<double>> DensityMatrix;
        RungeKutta<Operator<std::complex<double>>> RK_object;

    public:
        Simulation()
        {
            material = Material("/home/gcistaro/NEGF/tb_models/TBgraphene");

            auto& DensityMatrix = RK_object.get_Function();
            //setting up grids
            auto& MasterRgrid = DensityMatrix.get_Operator_R().get_MeshGrid();
            MasterRgrid = std::make_shared<MeshGrid<R>>(MeshGrid<R>(3.));//spherical grid for exponentially decaying DM
            
            auto& kgrid = DensityMatrix.get_Operator_k().get_MeshGrid();
            kgrid = std::make_shared<MeshGrid<k>>(std::move(fftPair<R, k>(*MasterRgrid)));//fftPair for k space
            
            auto Rfft = std::make_shared<MeshGrid<R>>(fftPair<k, R>(*kgrid));//Rspace as fourier space of kgrid
            print_grids();   
            material.print_info();

            //get U and Udagger from Hamiltonian
            SettingUp_EigenSystem();
            
            auto& Uk = Operator<std::complex<double>>::EigenVectors;
            auto& UkDagger = Operator<std::complex<double>>::EigenVectors_dagger;
            auto InitialCondition = [&](Operator<std::complex<double>>& DM)
            {
                DM.Operator_k.initialize(Uk.get_nblocks(), Uk.get_nrows(), Uk.get_ncols());
                DM.Operator_k.fill(0.);
                //filling matrix in Bloch gauge
                for(int ik=0; ik<Uk.get_nblocks(); ++ik){
                    for(int iband=0; iband<Uk.get_nrows(); iband++){
                        if(this->Band_energies[ik](iband) < FermiEnergy-threshold){
                            DM.Operator_k(ik, iband, iband) = 1.;
                        }
                    }
                }
                DM.lock_gauge(bloch);
                DM.lock_space(k);
                DM.go_to_wannier();
                DM.go_to_R();
            };
            
        };

        void SettingUp_EigenSystem()
        {
            //diagonalizing Hk on the k vectors of the grid of the Density matrix.
            auto& DensityMatrix = RK_object.get_Function();
            auto& MasterkGrid = DensityMatrix.get_Operator_k().get_MeshGrid()->get_mesh();
            std::cout << "MasterkGrid::: \n";
            for(int ik=0; ik<MasterkGrid.size(); ik++){
                std::cout << ik << MasterkGrid[ik].get("LatticeVectors");
            }
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
            auto& DensityMatrix = RK_object.get_Function();
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
        
};

