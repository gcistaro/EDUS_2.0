#include "BlockMatrix.hpp"


int main()
{
    mdarray<double,3> a({4,4,4});
    BlockMatrix<double> output;
    output.initialize(a);
    //std::cout << output[0];

    a=mdarray<double,3>({4,4,4});
    BlockMatrix<double> input1;
    input1.initialize(a);
    //std::cout << input1[0];
    
    a=mdarray<double,3>({4,4,4});
    BlockMatrix<double> input2;
    input2.initialize(a);
    //std::cout << input2[0];

    multiply(output, 1., input1, input2);

}