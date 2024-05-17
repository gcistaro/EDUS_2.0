std::function<void(BlockMatrix<std::complex<double>, R>&, double const&, BlockMatrix<std::complex<double>, R> const&)> SourceTerm = 
[&](BlockMatrix<std::complex<double>, R>& Output, const double& time, const BlockMatrix<std::complex<double>, R>& Input)
{
    PROFILE("RK::SourceTerm");
    //Von Neumann Equations: idp/dt = [H, p]

    Calculate_TDHamiltonian(time);
    Output.fill(0.*im);
    commutator(Output, -im, H.get_Operator_R(), Input);
    
    //R operator 
    #pragma omp parallel for schedule(dynamic)
    for(int iR=0; iR<Output.get_nblocks(); ++iR){
        auto& R = (*(Output.get_MeshGrid()))[iR].get("Cartesian");
        auto las = laser(time).get("Cartesian");
        auto DotProduct = las[0]*R[0] + las[1]*R[1] + las[2]*R[2];

        for(int irow=0; irow<Output.get_nrows(); ++irow){
            for(int icol=0; icol<Output.get_ncols(); ++icol){
                Output(iR,irow, icol) += im*Input(iR,irow, icol)*DotProduct;
            }
        }
    }
};


