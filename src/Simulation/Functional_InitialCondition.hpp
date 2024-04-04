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


