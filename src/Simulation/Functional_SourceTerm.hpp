std::function<void(Operator<std::complex<double>>&, double const&, Operator<std::complex<double>> const&)> 
SourceTerm = 
[&](Operator<std::complex<double>>& Output__, const double& time__, const Operator<std::complex<double>>& Input__)
{
    Output__.go_to_R(false);
    const_cast<Operator<std::complex<double>>&>(Input__).go_to_R(true);
    H_.go_to_R(false);
    H_.get_Operator_R().fill(0);
    
    //Output +=   (E.Nabla) * Input
    kgradient_.Calculate(1.+0.*im, Output__.get_Operator(Space::R), 
                        Input__.get_Operator(Space::R), 
                        setoflaser_(time__), true);
                        
    coulomb_.EffectiveHamiltonian( H_, Input__, true); 
    Output__.go_to_k(true);
    const_cast<Operator<std::complex<double>>&>(Input__).go_to_k(false);
    H_.go_to_k(true);


    // Output = -i * [ H, Input ]
    Calculate_TDHamiltonian(time__, false);
    auto& Output = Output__.get_Operator(Space::k);
    auto& Input = Input__.get_Operator(Space::k);
    auto& H = H_.get_Operator(Space::k);
    
    commutator(Output, -im, H, Input, false);
    //assert(Output_.is_hermitian());
};