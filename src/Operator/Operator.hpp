#ifndef OPERATOR_HPP
#define OPERATOR_HPP

#include "cassert"
#include "Operator/BlockMatrix.hpp"
#include "fftPair/fftPair.hpp"

class BandIndex
{
    private:
        size_t NumberOfBands;
        size_t oneDNumberOfBands;

        std::vector< std::pair<size_t,size_t> > RowIndexBoundary;

/*
    numbering:
    [0   1   2    3 ]
    [.   4   5    6 ]
    [.   .   7    8 ]
    [.   .   .    9 ]

    ......
    Total number of elements before row m:
    b + (b-1) + ... + (b-m-1) = \sum_{i=0}^{m-1} (b-i) = m*b-(m-1)*(m)/2  .... VALID FOR m>1
    This is the starting index of row m.

    The final index is
    StartingIndex + #elements in row = m*b-(m-1)*m/2 + (b-m-1) 

*/ 
    public:
    size_t StartingIndex(const auto& row_)
    {
        size_t StartingIndex = (row_ == 0) ? 0 
                                           :row_*NumberOfBands - row_*(row_-1)/2;
        return StartingIndex;
    };

    size_t RowIndex(const size_t& oneDband)
    {
        auto RowIndexIterator_ = std::find_if(RowIndexBoundary.begin(), RowIndexBoundary.end(),
                                             [&](const auto& RowPair){return oneDband >= RowPair.first && 
                                                                             oneDband <= RowPair.second;
                                                                      });
        return RowIndexIterator_ - RowIndexBoundary.begin();
    };

    BandIndex(){};
    BandIndex(const size_t& NumberOfBands__)
    {
        initialize(NumberOfBands__);
    };

    void initialize(const size_t& NumberOfBands__)
    {
        NumberOfBands = NumberOfBands__;
        oneDNumberOfBands = StartingIndex(NumberOfBands);
        RowIndexBoundary.resize(NumberOfBands);
        for(int i=0; i<NumberOfBands; ++i){
            RowIndexBoundary[i].first = StartingIndex(i);
            RowIndexBoundary[i].second = StartingIndex(i) + (NumberOfBands - i) - 1;
            //std::cout << "row " << i << " boundaries: [" << RowIndexBoundary[i].first << " , " << RowIndexBoundary[i].second << "]" << std::endl;
        }

    }
    size_t oneDband(const size_t& bnd1, const size_t& bnd2)
    {
        assert(bnd1 <= bnd2);
        assert(bnd1 <= NumberOfBands && bnd2 <= NumberOfBands);

        return RowIndexBoundary[bnd1].first + (bnd2-bnd1);
    };

    std::pair<size_t, size_t> twoDband(const size_t& oneDband__)
    {
        assert(oneDband__ < oneDNumberOfBands);
        std::pair<size_t, size_t> twoDband_;
        twoDband_.first = RowIndex(oneDband__);
        twoDband_.second = oneDband__ - RowIndexBoundary[twoDband_.first].first + twoDband_.first;
        return twoDband_;
    };

    size_t& get_oneDNumberOfBands(){ return oneDNumberOfBands; };
};




template < typename T=std::complex<double> >
class Operator
{
    private:
	    BlockMatrix<T,k> Operator_k;
        BlockMatrix<T,R> Operator_R;

        mdarray<std::complex<double>, 2> FTfriendly_Operator_R;
        mdarray<std::complex<double>, 2> FTfriendly_Operator_k;
        std::shared_ptr<MeshGrid<k>> FT_meshgrid_k;
        std::shared_ptr<MeshGrid<R>> FT_meshgrid_R;

        static BandIndex bandindex;
        BandGauge bandgauge;
        Space space;

        bool locked_bandgauge = false;
        bool locked_space = false;
        bool initialized_dft = false;
        bool initialized_fft = false;
        FourierTransform ft_;
        

    public:
        static BlockMatrix<T, k> temp_k;
        friend class Simulation;
        static std::shared_ptr<MeshGrid<R>> MeshGrid_Null;
        static BlockMatrix<T,k> EigenVectors;
        static BlockMatrix<T,k> EigenVectors_dagger;
        Operator() : Operator_k(BlockMatrix<T,k>()), Operator_R(BlockMatrix<T,R>()){};


        Operator(const Operator<T>& Op_) = default;
        Operator<T>& operator=(const Operator<T>& Op_) = default;
        
        Operator(Operator<T>&& Op_)  = default;
        Operator<T>& operator=(Operator<T>&& Op_) = default;

        BlockMatrix<T,k>& get_Operator_k()
        {
            return const_cast<BlockMatrix<T,k>&>(static_cast<const Operator<T>&>(*this).get_Operator_k());
        };       
        
        const BlockMatrix<T,k>& get_Operator_k() const
        {
            return Operator_k;
        };      
        
        BlockMatrix<T,R>& get_Operator_R()
        {
            return const_cast<BlockMatrix<T,R>&>(static_cast<const Operator<T>&>(*this).get_Operator_R());
        };      
        
        const BlockMatrix<T,R>& get_Operator_R() const
        {
            return Operator_R;
        };


        template<Space space>
        void initialize_fft(const MeshGrid<space>& MG, const size_t& nbnd)
        {
            //for now this is the only case implemented. it will be more general.
            //we enter in this function only once
            if(initialized_fft){
                return;
            }
            initialized_fft = true;
            
            Operator_R = BlockMatrix<std::complex<double>,R>(MG.get_mesh().size(), nbnd, nbnd);

            switch(space)
            {
                case R:
                {
                    auto& Rgrid = Operator_R.get_MeshGrid();
                    Rgrid = std::make_shared<MeshGrid<R>>(MG);
                    FT_meshgrid_k = std::make_shared<MeshGrid<k>>(fftPair<R,k>(MG));
                    FT_meshgrid_R = std::make_shared<MeshGrid<R>>(fftPair<k,R>(*FT_meshgrid_k));
                    break;
                }
                case k:
                {
                    //FT_meshgrid_R = std::make_shared<MeshGrid<R>>(fftPair<R,k>(MG));
                    //FT_meshgrid_k = std::make_shared<MeshGrid<k>>(fftPair<R,k>(*FT_meshgrid_R));
                    std::cout <<"Initialize fft from k Not yet implemented!\n";
                    exit(1);
                }
            }
            //mdarray<double,2> bare_mg({1,3});
            //bare_mg.fill(0);
            //MeshGrid_Null = std::make_shared<MeshGrid<R>>(bare_mg, "Cartesian");
            std::cout << "Calculate convolution index MG, FT_meshgrid_R, MeshGrid_null\n";
            MeshGrid<R>::Calculate_ConvolutionIndex1(MG , *FT_meshgrid_R, *MeshGrid_Null);
            std::cout << "Calculate convolution index MeshGrid_null, MG, MG\n";
            MeshGrid<R>::Calculate_ConvolutionIndex1(*MeshGrid_Null, MG, MG);

            //prove that is working
            std::ofstream MG1, MG2;
            MG1.open("Mg1.txt");
            MG2.open("Mg2.txt");
            for(int iR=0; iR<FT_meshgrid_R->get_mesh().size(); iR++){
                MG1 <<  (*FT_meshgrid_R)[iR].get("LatticeVectors");
            }
            for(int iR=0; iR<MG.get_mesh().size(); iR++){
                auto& ci = MeshGrid<R>::ConvolutionIndex1[{MG.get_id(), FT_meshgrid_R->get_id(), MeshGrid_Null->get_id()}];
                if(ci(iR, 0) != -1){
                    MG2 << (*FT_meshgrid_R)[ci(iR,0)].get("LatticeVectors"); 
                }
            }
            MG1.close();
            MG2.close();
            //endofproof
            
            auto nk = FT_meshgrid_R->get_mesh().size();

            Operator_k = BlockMatrix<std::complex<double>,k>(nk, nbnd, nbnd);
            auto& kgrid = Operator_k.get_MeshGrid();
            kgrid = FT_meshgrid_k;
            bandindex.initialize(nbnd);
            

            
            FTfriendly_Operator_k = mdarray<std::complex<double>, 2>({nbnd*(nbnd+1)/2, nk});
            FTfriendly_Operator_R = mdarray<std::complex<double>, 2>({nbnd*(nbnd+1)/2, nk});

            //use convolution index for shuffle index.
            std::vector<int> Dimensions(3);
            for(int ix=0; ix<3; ix++){
                Dimensions[ix] = FT_meshgrid_k->get_Size()[ix];
            }

            ft_.initialize(FTfriendly_Operator_k, FTfriendly_Operator_R, Dimensions);
            //shuffle_to_fft();
        };


        void dft(const std::vector<Coordinate<k>>& path, const int& sign)
        {
            //for(auto& k : path){
            //    auto& kk = k.get("LatticeVectors");
            //    std::cout << kk[0] << " " << kk[1] << " " << kk[2] << std::endl;
            //}
            initialize_dft();
            execute_dft(path, sign);
            shuffle_to_RK();
        }

/*
 * "initialize_dft()"
 * The following function is used to reshuffle the operator in such a way that R is the first index
 * and can be easily used for doing a fft. In the same time, we put as R mesh of the dft the one of 
 * Operator_R.
*/
        void initialize_dft()
        {
            //we enter in this function only once
            if(initialized_dft){
                return;
            }
            initialized_dft = true;
            auto nbnd = Operator_R.get_nrows();
            bandindex.initialize(nbnd);
            //initialize input of fft
            FTfriendly_Operator_R = mdarray<std::complex<double>, 2>({nbnd*(nbnd+1)/2, Operator_R.get_nblocks()});
            for(int iR=0; iR<Operator_R.get_nblocks(); iR++){
                for(int ibnd1=0; ibnd1<nbnd; ++ibnd1){
                    for(int ibnd2=ibnd1; ibnd2<nbnd; ++ibnd2){
                        FTfriendly_Operator_R(static_cast<int>(bandindex.oneDband(ibnd1, ibnd2)), iR) = 
                                                                            Operator_R(iR, ibnd1, ibnd2);
                    }
                }
            }
            std::vector< std::vector<double> > Mesh_FT;
            auto& mesh_operator = Operator_R.get_MeshGrid()->get_mesh();
            Mesh_FT.resize(mesh_operator.size());        

            for(int im=0; im<mesh_operator.size(); im++){
                Mesh_FT[im].resize(3);
                Mesh_FT[im][0] = mesh_operator[im].get("LatticeVectors")[0];
                Mesh_FT[im][1] = mesh_operator[im].get("LatticeVectors")[1];
                Mesh_FT[im][2] = mesh_operator[im].get("LatticeVectors")[2];
                //std::cout <<"MESH_FT: "<< Mesh_FT[im][0] <<" " << Mesh_FT[im][1] << " " << Mesh_FT[im][2] << std::endl;
                //std::cout << Mesh_FT[im][0] << " " << Mesh_FT[im][1] << Mesh_FT[im][2] << std::endl;
            }
            ft_.initialize(FTfriendly_Operator_R, Mesh_FT);
            std::cout << "ft is initialized.\n";
        }


        void execute_dft(const std::vector<Coordinate<k>>& path, const int& sign)
        {
            std::vector< std::vector<double>> path_bare(path.size());
            for(int i=0; i<path.size(); ++i){
                path_bare[i].resize(3);
                path_bare[i][0] = path[i].get("LatticeVectors")[0];
                path_bare[i][1] = path[i].get("LatticeVectors")[1];
                path_bare[i][2] = path[i].get("LatticeVectors")[2];
            }
            FTfriendly_Operator_k = ft_.dft(path_bare, +1);
        }


        void shuffle_to_RK()
        {
            Operator_k.initialize(FTfriendly_Operator_k.get_Size(1), Operator_R.get_nrows(), Operator_R.get_ncols());
            for(int ik=0; ik<Operator_k.get_nblocks(); ++ik){
                for(int ibnd1=0; ibnd1<Operator_k.get_nrows(); ++ibnd1){
                    for(int ibnd2=ibnd1; ibnd2<Operator_k.get_nrows(); ++ibnd2){ 
                        Operator_k(ik, ibnd1, ibnd2) = FTfriendly_Operator_k(static_cast<int>(bandindex.oneDband(ibnd1,ibnd2)),ik);            
                        Operator_k(ik, ibnd2, ibnd1) = std::conj(Operator_k(ik,ibnd1, ibnd2));
        		    }
                }
            }
        }

        void shuffle_to_fft_R()
        {
            auto ci = MeshGrid<R>::get_ConvolutionIndex1(*Operator_R.get_MeshGrid() , *FT_meshgrid_R, *MeshGrid_Null);
            for(int iR=0; iR<Operator_R.get_nblocks(); iR++){
                assert(ci(iR,0) != -1);
                for(int ibnd1=0; ibnd1<Operator_R.get_nrows(); ++ibnd1){
                    for(int ibnd2=ibnd1; ibnd2<Operator_R.get_ncols(); ++ibnd2){
                        //std::cout  << " R " << iR << "ibnd1 " << ibnd1  << "ibnd2 "<< ibnd2  << " ci(iR,0) "<<   ci(iR,0) << std::endl;
                        FTfriendly_Operator_R(static_cast<int>(bandindex.oneDband(ibnd1, ibnd2)), ci(iR,0)) = Operator_R(iR, ibnd1, ibnd2);
                    }
                }
            }
        }

        void shuffle_to_fft_k()
        {
            for(int ik=0; ik<Operator_k.get_nblocks(); ik++){
                for(int ibnd1=0; ibnd1<Operator_k.get_nrows(); ++ibnd1){
                    for(int ibnd2=ibnd1; ibnd2<Operator_k.get_ncols(); ++ibnd2){
                        FTfriendly_Operator_k(static_cast<int>(bandindex.oneDband(ibnd1, ibnd2)), ik) = Operator_k(ik, ibnd1, ibnd2);
                        //std::cout << bandindex.oneDband(ibnd1, ibnd2) << " " << ik;
                        //std::cout << Operator_k(ik, ibnd1, ibnd2) << std::endl;
                    }
                }
            }
        }

        void shuffle_from_fft_R()
        {
            auto ci = MeshGrid<R>::get_ConvolutionIndex1(*Operator_R.get_MeshGrid(), *FT_meshgrid_R, *MeshGrid_Null);
            auto ciminus = MeshGrid<R>::get_ConvolutionIndex1(*MeshGrid_Null, *Operator_R.get_MeshGrid(), *Operator_R.get_MeshGrid());

            for(int iR=0; iR<Operator_R.get_nblocks(); iR++){
                assert(ci(iR,0) != -1);
                for(int ibnd1=0; ibnd1<Operator_R.get_nrows(); ++ibnd1){
                    for(int ibnd2=ibnd1; ibnd2<Operator_R.get_ncols(); ++ibnd2){
                        //std::cout  << " R " << iR << "ibnd1 " << ibnd1  << "ibnd2 "<< ibnd2  << " ci(iR,0) "<<   ci(iR,0) << std::endl;
                        Operator_R(iR, ibnd1, ibnd2) = FTfriendly_Operator_R(static_cast<int>(bandindex.oneDband(ibnd1, ibnd2)), ci(iR,0));
                    }
                }
            }

            //this  part is not working! Check ciminus
            for(int iR=0; iR<Operator_R.get_nblocks(); iR++){
                assert(ciminus(0,iR) != -1);
                for(int ibnd1=0; ibnd1<Operator_R.get_nrows(); ++ibnd1){
                    for(int ibnd2=ibnd1+1; ibnd2<Operator_R.get_ncols(); ++ibnd2){
                        //std::cout << iR << " " << (*(Operator_R.get_MeshGrid()))[iR] <<
                        Operator_R(iR, ibnd2, ibnd1) = std::conj(Operator_R(ciminus(0,iR), ibnd1, ibnd2));
                    }
                }
            }
        }

        void shuffle_from_fft_k()
        {
            for(int ik=0; ik<Operator_k.get_nblocks(); ik++){
                for(int ibnd1=0; ibnd1<Operator_k.get_nrows(); ++ibnd1){
                    for(int ibnd2=ibnd1; ibnd2<Operator_k.get_ncols(); ++ibnd2){
                        Operator_k(ik, ibnd1, ibnd2) = FTfriendly_Operator_k(static_cast<int>(bandindex.oneDband(ibnd1, ibnd2)), ik);
                        Operator_k(ik, ibnd2, ibnd1) = conj(Operator_k(ik, ibnd1, ibnd2));
                    }
                }
            }
        }


        void go_to_wannier()
        {
            assert(locked_bandgauge);
            //assert(space == k);
            if (bandgauge == wannier){
                return;
            }
            auto temp_k = Operator_k;
            temp_k.fill(0); 

            multiply(temp_k, std::complex<double>(1.), EigenVectors, Operator_k);
            Operator_k.fill(0);
            multiply(Operator_k, std::complex<double>(1.), temp_k, EigenVectors_dagger);
            bandgauge = wannier;
        };

        void go_to_bloch()
        {
            assert(locked_bandgauge);
            //assert(space == k);
            if(bandgauge == bloch){
                return;
            }
            auto temp_k = Operator_k;
            temp_k.fill(0); 
            multiply(temp_k, std::complex<double>(1.), EigenVectors_dagger, Operator_k);
            Operator_k.fill(0);
            multiply(Operator_k, std::complex<double>(1.), temp_k, EigenVectors);
            bandgauge = bloch;
        };

        void go_to_R()
        {
            assert(locked_space && initialized_fft);
            //if(space == R){
            //    return;
            //}
            shuffle_to_fft_k();
            ft_.fft(-1);       
            shuffle_from_fft_R();  

            //space = R;
        }

        void go_to_k()
        {
            assert(locked_space && initialized_fft);
            //if(space == k){
            //    return;
            //}
            shuffle_to_fft_R();
            ft_.fft(+1);          
            shuffle_from_fft_k();  

            //space = k;
        }

        void lock_gauge(const BandGauge& bandgauge__)
        {
            bandgauge = bandgauge__;
            locked_bandgauge = true;
        };

        void lock_space(const Space& space__)
        {
            space = space__;
            locked_space = true;
        };

};





template < typename T>
BandIndex Operator<T>::bandindex;

template < typename T>
BlockMatrix<T,k> Operator<T>::EigenVectors;

template < typename T>
BlockMatrix<T,k> Operator<T>::EigenVectors_dagger;


template < typename T>
BlockMatrix<T,k> Operator<T>::temp_k;

template < typename T>
std::shared_ptr<MeshGrid<R>> Operator<T>::MeshGrid_Null = std::make_shared<MeshGrid<R>>(MeshGrid<R>(std::vector<Coordinate<R>>({Coordinate<R>(0,0,0)})));


#endif