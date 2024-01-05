#include "BlockMatrix.hpp"


template < typename T=std::complex<double>>
class Operator
{
    private:
	BlockMatrix<T,k> Operator_k;
        BlockMatrix<T,R> Operator_R;

        static BlockMatrix<T,k> EigenVectors;
        static BlockMatrix<T,k> EigenVectors_dagger;
        
        //fft parameters
        //std::shared_ptr<MeshGrid<R>> R_fft;
        //std::map< std::pair<int,int>, int> R_shuffle;

        enum BandGauge{Bloch, Wannier};

        

    public:
        Operator() : Operator_k(BlockMatrix<T,k>()), Operator_R(BlockMatrix<T,R>())
	{};//this is just to make them callable and copyable
       
        Operator(const Operator<T>& Op_) {*this = Op_;}
        Operator<T>& operator=(const Operator<T>& Op_)
	{
		std::cout << "Copying operator..\n";
	    Operator_k = Op_.Operator_k;
	    Operator_R = Op_.Operator_R;
            return *this;
	}

        Operator(Operator<T>&& Op_)  {*this = Op_;}
        Operator<T>& operator=(Operator<T>&& Op_)
	{
		std::cout << "Moving operator...\n";
	    Operator_k = Op_.Operator_k;
	    Operator_R = Op_.Operator_R;
            return *this;
	}
        	
        BlockMatrix<T,k>& get_Operator_k()
        {
            return const_cast<BlockMatrix<T,k>&>(static_cast<const Operator<T>&>(*this).get_Operator_k());
        }

        const BlockMatrix<T,k>& get_Operator_k() const
        {
            return Operator_k;
        }


        BlockMatrix<T,R>& get_Operator_R()
        {
            return const_cast<BlockMatrix<T,R>&>(static_cast<const Operator<T>&>(*this).get_Operator_R());
        }

        const BlockMatrix<T,R>& get_Operator_R() const
        {
            return Operator_R;
        }


};
