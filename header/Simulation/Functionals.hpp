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

std::function<void(BlockMatrix<std::complex<double>, k>&)> InitialCondition_k = [&](BlockMatrix<std::complex<double>, k>& DM)
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
    //DensityMatrix.go_to_R();
};

std::function<void(BlockMatrix<std::complex<double>, R>&, double const&, BlockMatrix<std::complex<double>, R> const&)> SourceTerm = 
[&](BlockMatrix<std::complex<double>, R>& Output, const double& time, const BlockMatrix<std::complex<double>, R>& Input)
{
    PROFILE("RK::SourceTerm");
    //Von Neumann Equations: idp/dt = [H, p]
    Output.fill(0.*im);
    Calculate_TDHamiltonian(time);
    convolution1(Output, -im, H.get_Operator_R(), Input);
    convolution1(Output, +im, Input, H.get_Operator_R() );
    
    std::cout << "Max of Input: " << max(Input) << std::endl;
    std::cout << "Max of Output: " << max(Output) << std::endl;
    auto& Rmesh = Output.get_MeshGrid()->get_mesh();
    for(auto& R : Rmesh){
        std::cout << R.get("LatticeVectors")[0] << " " <<R.get("LatticeVectors")[1] << " " << R.get("LatticeVectors")[2]<< std::endl;
    }
};


std::function<void(BlockMatrix<std::complex<double>, k>&, double const&, BlockMatrix<std::complex<double>, k> const&)> SourceTerm_k = 
[&](BlockMatrix<std::complex<double>, k>& Output, const double& time, const BlockMatrix<std::complex<double>, k>& Input)
{
    PROFILE("RK::SourceTerm");
    //Von Neumann Equations: idp/dt = [H, p]
    Output.fill(0.*im);
    Calculate_TDHamiltonian(time);
    multiply(Output, -im, material.H.get_Operator_k(), Input);
    //material.H.go_to_bloch();
    //std::cout << "materialH: " << material.H.get_Operator_k();
    
    multiply(Output, +im, Input, material.H.get_Operator_k(), 1.+0.*im );
    
    std::cout << "Max of Input: " << max(Input) << std::endl;
    std::cout << "Max of Output: " << max(Output) << std::endl;
    std::cout << "materialH: " << material.H.get_Operator_k();
    std::cout << "Input: " << Input;
    std::cout << "Output: " << Output;
};
