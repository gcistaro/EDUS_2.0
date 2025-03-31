std::function<void(Operator<std::complex<double>>&, double const&, Operator<std::complex<double>> const&)> 
SourceTerm = 
[&](Operator<std::complex<double>>& Output__, const double& time__, const Operator<std::complex<double>>& Input__)
{
    
    /* Gradient term:    Output +=   (E.Nabla) * Input */
    Output__.go_to_R(false);
    const_cast<Operator<std::complex<double>>&>(Input__).go_to_R(true);
    kgradient_.Calculate(1.+0.*im, Output__.get_Operator(SpaceOfPropagation_Gradient_), 
                        Input__.get_Operator(SpaceOfPropagation_Gradient_), 
                        setoflaser_(time__), true);

    /* Coulomb interaction */
    H_.go_to_R(false);
    H_.get_Operator_R().fill(0);
    coulomb_.EffectiveHamiltonian( H_, Input__, true); 

    if ( SpaceOfPropagation_Gradient_ == R ) {
        Output__.go_to_k(true);
        const_cast<Operator<std::complex<double>>&>(Input__).go_to_k(false);
    }
    H_.go_to_k(true);


    // Output = -i * [ H, Input ]
    Calculate_TDHamiltonian(time__, false);
    auto& Output = Output__.get_Operator(Space::k);
    auto& Input = Input__.get_Operator(Space::k);
    auto& H = H_.get_Operator(Space::k);
    commutator(Output, -im, H, Input, false);

};