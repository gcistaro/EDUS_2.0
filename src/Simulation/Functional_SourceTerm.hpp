std::function<void(Operator<std::complex<double>>&, double const&, Operator<std::complex<double>> const&)> 
SourceTerm = 
[&](Operator<std::complex<double>>& Output, const double& time, const Operator<std::complex<double>>& Input)
{
    Output.go_to_R(false);
    const_cast<Operator<std::complex<double>>&>(Input).go_to_R(true);
    H.go_to_R(false);
    H.get_Operator_R().fill(0);
    
    //Output +=   (E.Nabla) * Input
    kgradient.Calculate(1.+0.*im, Output.get_Operator(Space::R), 
                        Input.get_Operator(Space::R), 
                        setoflaser(time), true);
                        
    coulomb.EffectiveHamiltonian( H, Input, true); 
    Output.go_to_k(true);
    const_cast<Operator<std::complex<double>>&>(Input).go_to_k(false);
    H.go_to_k(true);

    //write Hamiltonian
    if( PrintObservables(time) ) 
        H.get_Operator_k().write_h5("SelfEnergy", std::to_string(get_it(time)) );

    // Output = -i * [ H, Input ]
    Calculate_TDHamiltonian(time, false);
    auto& Output_ = Output.get_Operator(Space::k);
    auto& Input_ = Input.get_Operator(Space::k);
    auto& H_ = H.get_Operator(Space::k);
    
    commutator(Output_, -im, H_, Input_, false);
    //assert(Output_.is_hermitian());
};