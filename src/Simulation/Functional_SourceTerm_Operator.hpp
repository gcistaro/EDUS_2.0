std::function<void(Operator<std::complex<double>>&, double const&, Operator<std::complex<double>> const&)> 
SourceTerm = 
[&](Operator<std::complex<double>>& Output, const double& time, const Operator<std::complex<double>>& Input)
{
    // dP/dt = -i( [H0+E.r, P] - (E.R)P )
    // dP/dt <- Output  
    // P <- Input
    PROFILE("RK::SourceTerm");
    auto& Input_  = Input.get_Operator(SpaceOfPropagation);
    auto& Output_ = Output.get_Operator(SpaceOfPropagation);

    // H = H0 + E.r
    Calculate_TDHamiltonian(time);

    // Output = -i * [ H, Input ]
    Output_.fill(0.*im);
    commutator(Output_, -im, H.get_Operator(SpaceOfPropagation), Input_);
    
    // Output += i * (E.R) * Input      (in R space)
    if ( SpaceOfPropagation == R ) {
        #pragma omp parallel for schedule(dynamic)
        for(int iR=0; iR<Output_.get_nblocks(); ++iR) {        
            auto& R = (*(Output_.get_MeshGrid()))[iR];
            auto imDotProduct = im*R.dot( laser(time) );

            Output_[iR] += imDotProduct * Input_[iR];
        }
    } else if ( SpaceOfPropagation == k ) {
        kgradient.Calculate(Output_, Input_, laser(time), false);
    }
};