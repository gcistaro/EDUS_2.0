std::function<void(Operator<std::complex<double>>&, double const&, Operator<std::complex<double>> const&)> 
SourceTerm = 
[&](Operator<std::complex<double>>& Output, const double& time, const Operator<std::complex<double>>& Input)
{
    // dP/dt = -i( [H0+E.r, P] - (E.R)P )
    // dP/dt <- Output  
    // P <- Input
    PROFILE("RK::SourceTerm");

    H.go_to_R(false);
    Output.go_to_R();
    const_cast<Operator<std::complex<double>>&>(Input).go_to_R();

    // Output +=   (E.Nabla) * Input
    PROFILE_START("i*(E.R)*Input");
    kgradient.Calculate(Output.get_Operator(Space::R), 
                        Input.get_Operator(Space::R), 
                        laser(time), false);
    PROFILE_STOP("i*(E.R)*Input");

    PROFILE_START("Coulomb");
    coulomb.EffectiveHamiltonian(H, Input, true);
    PROFILE_STOP("Coulomb");

    Output.go_to_k();
    const_cast<Operator<std::complex<double>>&>(Input).go_to_k(false);
    H.go_to_k();

    // H = H0 + E.r
    Calculate_TDHamiltonian(time, !coulomb.get_DoCoulomb() );
    
    // Output = -i * [ H, Input ]
    auto& Output_ = Output.get_Operator(Space::k);
    auto& Input_ = Input.get_Operator(Space::k);
    auto& H_ = H.get_Operator(Space::k);
    
    Output_.fill(0.*im);
    commutator(Output_, -im, H_, Input_);
    //assert(Output_.is_hermitian());
    //std::cout << *max(Output_) << std::endl;
};