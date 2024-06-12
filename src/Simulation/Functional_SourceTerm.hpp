std::function<void(Operator<std::complex<double>>&, double const&, Operator<std::complex<double>> const&)> 
SourceTerm = 
[&](Operator<std::complex<double>>& Output, const double& time, const Operator<std::complex<double>>& Input)
{
    // dP/dt = -i( [H0+E.r, P] - (E.R)P )
    // dP/dt <- Output  
    // P <- Input
    PROFILE("RK::SourceTerm");

    // H = H0 + E.r
    Calculate_TDHamiltonian(time);

    // Output = -i * [ H, Input ]
    auto& Output_ = Output.get_Operator(SpaceOfPropagation);
    auto& Input_ = Input.get_Operator(SpaceOfPropagation);
    auto& H_ = H.get_Operator(SpaceOfPropagation);
    
    Output_.fill(0.*im);
    commutator(Output_, -im, H_, Input_);
    //assert(Output_.is_hermitian());


    if( SpaceOfPropagation == R ) {
    // Output += i * (E.R) * Input
        #pragma omp parallel for schedule(dynamic)
        for(int iR=0; iR<Output_.get_nblocks(); ++iR){        
            auto& R = (*(Output_.get_MeshGrid()))[iR];
            auto imDotProduct = im*R.dot( laser(time) );
            Output_[iR] += imDotProduct * Input_[iR];
        }
    }
    else if ( SpaceOfPropagation == k ) {
    // Output +=   (E.Nabla) * Input
        kgradient.Calculate(Output_, Input_, 
                            laser(time), false);
    }
    //assert(Output_.is_hermitian());

    //std::cout << *max(Output_) << std::endl;
};