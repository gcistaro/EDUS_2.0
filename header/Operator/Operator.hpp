#include "BlockMatrix.hpp"


template < typename T=std::complex<double>, Space space >
class Operator : public BlockMatrix<T, space>
{
    private:
        BlockMatrix<T,k> Operator_k;
        BlockMatrix<T,R> Operator_R;

        static BlockMatrix<T,k> EigenVectors;
        static BlockMatrix<T,k> EigenVectors_dagger;



    public:
        Operator(){};
        Operator(const Space space_) : space(space_){};
        
        void Transform(Type& type_)
        {
            if(type == type_){
                return;
            }
            
            //all the following must be replaced with spFFT
            if(type == k){
                sign = +1;
            }
            else{
                sign = -1;
            }

        }
};