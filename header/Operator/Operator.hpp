#include "BlockMatrix.hpp"


template <typename T>
class Operator : public BlockMatrix<T>
{
    private:
        static BlockMatrix<T> EigenVectors;
        static BlockMatrix<T> EigenVectors_dagger;

        Type type;
    
    public:
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