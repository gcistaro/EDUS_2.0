std::function<void(Operator<std::complex<double>>&)> 
InitialCondition = 
[&](Operator<std::complex<double>>& DM)
{
    PROFILE("RK::InitialCondition");
    DM.get_Operator_k().fill(0.);
    
    //filling matrix in Bloch-k gauge
    for(int ik=0; ik<Uk.get_nblocks(); ++ik){
        for(int iband=0; iband<Uk.get_nrows(); iband++){
            if(this->Band_energies[ik](iband) < FermiEnergy-threshold){
                DM.get_Operator_k()(ik, iband, iband) = 1.;
            }
        }
    }
    DM.lock_gauge(bloch);
    DM.lock_space(k);
    DM.go_to_wannier();
    DM.go_to_R();
};


