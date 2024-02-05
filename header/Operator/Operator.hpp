#include "cassert"
#include "BlockMatrix.hpp"
#include "fftPair.hpp"

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
            std::cout << "row " << i << " boundaries: [" << RowIndexBoundary[i].first << " , " << RowIndexBoundary[i].second << "]" << std::endl;
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

        static BlockMatrix<T,k> EigenVectors;
        static BlockMatrix<T,k> EigenVectors_dagger;
        static BandIndex bandindex;
        enum BandGauge{Bloch, Wannier};

        

    public:
        Operator() : Operator_k(BlockMatrix<T,k>()), Operator_R(BlockMatrix<T,R>()){};

        Operator(const Operator<T>& Op_) {*this = Op_;}
        
        Operator<T>& operator=(const Operator<T>& Op_)
        {
            Operator_k = Op_.Operator_k;
            Operator_R = Op_.Operator_R;
            return *this;
        };
        
        Operator(Operator<T>&& Op_)  {*this = Op_;}
        Operator<T>& operator=(Operator<T>&& Op_)
        {
        	std::cout << "Moving operator...\n";
            Operator_k = Op_.Operator_k;
            Operator_R = Op_.Operator_R;
                return *this;
        };

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


        void initialize_fft()
        {
            auto nbnd = Operator_R.get_nrows();
            bandindex.initialize(nbnd);
            
            FT_meshgrid_k = std::make_shared<MeshGrid<k>>(fftPair<R,k>(*Operator_R.get_MeshGrid()));
            FT_meshgrid_R = std::make_shared<MeshGrid<R>>(fftPair<k,R>(*FT_meshgrid_k));

            auto nk = FT_meshgrid_R->get_mesh().size();
            
            FTfriendly_Operator_k = mdarray<std::complex<double>, 2>({nbnd*(nbnd+1)/2, nk});
            FTfriendly_Operator_R = mdarray<std::complex<double>, 2>({nbnd*(nbnd+1)/2, nk});

            //use convolution index for shuffle index.
            mdarray<double,2> bare_mg({1,3});
            bare_mg.fill(0);

            auto MeshGrid_Null = std::make_shared<MeshGrid<R>>(bare_mg, "Cartesian");
            MeshGrid<R>::Calculate_ConvolutionIndex(*Operator_R.get_MeshGrid() , *FT_meshgrid_R, *MeshGrid_Null);
            auto ci = MeshGrid<R>::get_ConvolutionIndex(*Operator_R.get_MeshGrid() , *FT_meshgrid_R, *MeshGrid_Null);
            for(int iR=0; iR<Operator_R.get_nblocks(); iR++){
                for(int ibnd1=0; ibnd1<nbnd; ++ibnd1){
                    for(int ibnd2=ibnd1; ibnd2<nbnd; ++ibnd2){
                        std::cout  << " R " << iR << "ibnd1 " << ibnd1  << "ibnd2 "<< ibnd2  << " ci(iR,0) "<<   ci(iR,0) << std::endl;
                        FTfriendly_Operator_R(bandindex.oneDband(ibnd1, ibnd2), ci(iR,0)) = Operator_R(iR, ibnd1, ibnd2);
                    }
                }
            }
        };
};

template < typename T>
BandIndex Operator<T>::bandindex;