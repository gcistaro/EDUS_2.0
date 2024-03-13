std::function<void(BlockMatrix<std::complex<double>, R>&)> InitialCondition = [&](BlockMatrix<std::complex<double>, R>& DM)
{
    PROFILE("RK::InitialCondition");
    //DM.Operator_k.initialize(Uk.get_nblocks(), Uk.get_nrows(), Uk.get_ncols());
    DensityMatrix.Operator_k.fill(0.);
    //filling matrix in Bloch gauge
    for(int ik=0; ik<Uk.get_nblocks(); ++ik){
        for(int iband=0; iband<Uk.get_nrows(); iband++){
            if(this->Band_energies[ik](iband) < FermiEnergy-threshold){
                DensityMatrix.Operator_k(ik, iband, iband) = 1.;
            }
        }
    }
    DensityMatrix.lock_gauge(bloch);
    DensityMatrix.lock_space(k);
    DensityMatrix.go_to_wannier();
    DensityMatrix.go_to_R();


};


std::function<void(BlockMatrix<std::complex<double>, R>&, double const&, BlockMatrix<std::complex<double>, R> const&)> SourceTerm = 
[&](BlockMatrix<std::complex<double>, R>& Output, const double& time, const BlockMatrix<std::complex<double>, R>& Input)
{
    PROFILE("RK::SourceTerm");
    //Von Neumann Equations: idp/dt = [H, p]

    Calculate_TDHamiltonian(time);
    Output.fill(0.*im);
    commutator(Output, -im, H.get_Operator_R(), Input);
    
    //R operator 
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

    //std::cout << "Max of Input: " << *max(Input) << std::endl;
    //std::cout << "H[0]: " << H.get_Operator_R()[0] << std::endl;
    //std::cout << "H[1]: " << H.get_Operator_R()[1] << std::endl;
    //std::cout << "Max of Output: " << *(max(Output))  << std::endl;

};


