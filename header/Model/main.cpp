#include "Model.hpp"

int main()
{
    //Model model("Convolution");
    Material model("TBgraphene");

    auto& H = model.H;
    H.initialize_fft();
    //auto& r0R = model.r[0].get_Operator_R();
    //auto& r1R = model.r[1].get_Operator_R();
    //auto& r2R = model.r[2].get_Operator_R();
    //mdarray<std::complex<double>,3> bare_output({HR.get_nblocks(), HR.get_nrows(), HR.get_ncols()});
    //std::cout << bare_output.get_Size(0) << bare_output.get_Size(1) << bare_output.get_Size(2) << std::endl; 
    //BlockMatrix<std::complex<double>,R> Output;
    //Output.initialize(bare_output);
    //std::cout << Output.get_nblocks() << Output.get_nrows() << Output.get_ncols() << std::endl; 
    //Output.set_MeshGrid((*r2R.get_MeshGrid()));
    //convolution(Output, 1., r0R, r1R);
//
    //std::ofstream fp;
    //auto mg = *r2R.get_MeshGrid();
    //fp.open("Output.txt");
//
    //for(int iR=0; iR<mg.get_TotalSize(); ++iR){
    //    fp << Output[iR](0,0).real() << std::endl;
    //}
}
