std::function<void(BlockMatrix<std::complex<double>, R>&, double const&, BlockMatrix<std::complex<double>, R> const&)> SourceTerm = 
[&](BlockMatrix<std::complex<double>, R>& Output, const double& time, const BlockMatrix<std::complex<double>, R>& Input)
{
    PROFILE("RK::SourceTerm");
    //Von Neumann Equations: idp/dt = [H, p]

    Calculate_TDHamiltonian(time);
    Output.fill(0.*im);
    commutator(Output, -im, H.get_Operator_R(), Input);
    //convolution1(Output, -im, H.get_Operator_R(), Input);
    //convolution1(Output, +im, Input, H.get_Operator_R() );
    
    //R operator 
    for(int iR=0; iR<Output.get_nblocks(); ++iR){
        auto& R = (*(Output.get_MeshGrid()))[iR].get("Cartesian");
        auto las = laser(time).get("Cartesian");
        auto DotProduct = las[0]*R[0] + las[1]*R[1] + las[2]*R[2];
        //auto R = (*(Output.get_MeshGrid()))[iR].get("LatticeVectors");
        //auto R0 = R[0]-20./2.;
        //auto R1 = R[1]-20./2.;
        //auto R_ = Coordinate<k>(R0,R1,0,"LatticeVectors");
        //auto& RR = R_.get("Cartesian");
        //auto DotProduct = las[0]*RR[0] + las[1]*RR[1] + las[2]*RR[2];


        for(int irow=0; irow<Output.get_nrows(); ++irow){
            for(int icol=0; icol<Output.get_ncols(); ++icol){
                Output(iR,irow, icol) += -im*Input(iR,irow, icol)*DotProduct;
            }
        }
    }

    //std::cout << "Max of Input: " << *max(Input) << std::endl;
    //std::cout << "H[0]: " << H.get_Operator_R()[0] << std::endl;
    //std::cout << "H[1]: " << H.get_Operator_R()[1] << std::endl;
    //std::cout << "Max of Output: " << *(max(Output))  << std::endl;

};


