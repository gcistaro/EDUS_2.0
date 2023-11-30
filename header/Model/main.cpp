#include "Model.hpp"

int main()
{
    Model model("Convolution");

    mdarray<std::complex<double>,3> bare_output({model.H.get_nblocks(), model.H.get_nrows(), model.H.get_ncols()});
    std::cout << bare_output.get_Size(0) << bare_output.get_Size(1) << bare_output.get_Size(2) << std::endl; 
    BlockMatrix<std::complex<double>,R> Output;
    Output.initialize(bare_output);
    std::cout << Output.get_nblocks() << Output.get_nrows() << Output.get_ncols() << std::endl; 
    Output.set_MeshGrid((*model.r[2].get_MeshGrid()));
    convolution(Output, 1., model.r[0], model.r[1]);

    std::ofstream fp;
    auto mg = *model.r[2].get_MeshGrid();
    fp.open("Output.txt");

    for(int iR=0; iR<model.r[2].get_MeshGrid()->get_TotalSize(); ++iR){
        fp << Output[iR](0,0).real() << std::endl;
    }
}