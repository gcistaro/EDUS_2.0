std::function<void(Operator<std::complex<double>>&)> 
InitialCondition = 
[&](Operator<std::complex<double>>& DM__)
{
    PROFILE("RK::InitialCondition");
    DM__.get_Operator_k().fill(0.);
    
    //filling matrix in Bloch-k gauge
    for(int ik=0; ik<Uk.get_nblocks(); ++ik){
        for(int iband=0; iband<ctx_->cfg().filledbands(); iband++){
        //for(int iband=0; iband<Uk.get_nrows(); iband++){
        //    if(this->Band_energies[ik](iband) < FermiEnergy-threshold){
                DM__.get_Operator_k()(ik, iband, iband) = 1.;
        //    }
        }
    }
    DM__.lock_gauge(bloch);
    DM__.lock_space(k);
    DM__.go_to_wannier();

    //assert(DM.get_Operator_k().is_hermitian());
    if(SpaceOfPropagation_ == R) DM__.go_to_R();
};


