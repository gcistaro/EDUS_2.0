std::function<void(Operator<std::complex<double>>&, double const&, Operator<std::complex<double>> const&)> 
SourceTerm = 
[&](Operator<std::complex<double>>& Output, const double& time, const Operator<std::complex<double>>& Input)
{
    Output.go_to_R(false);
    const_cast<Operator<std::complex<double>>&>(Input).go_to_R(true);
    H.go_to_R(false);
    H.get_Operator_R().fill(0);
    
    //Output +=   (E.Nabla) * Input
    PROFILE_START("i*(E.R)*Input");
    kgradient.Calculate(1.+0.*im, Output.get_Operator(Space::R), 
                        Input.get_Operator(Space::R), 
                        laser(time), true);
    PROFILE_STOP("i*(E.R)*Input");

    PROFILE_START("Coulomb");
    coulomb.EffectiveHamiltonian( H, Input, true); 
    PROFILE_STOP("Coulomb");

//std::cout << "HRR" << H.get_Operator(R) << std::endl<< std::endl<< std::endl<< std::endl<< std::endl<< std::endl<< std::endl;
    Output.go_to_k(true);
    const_cast<Operator<std::complex<double>>&>(Input).go_to_k(false);
    H.go_to_k(true);
//std::cout << "HK: " << H.get_Operator(k) << std::endl<< std::endl<< std::endl<< std::endl<< std::endl<< std::endl<< std::endl;

    // Output = -i * [ H, Input ]
    Calculate_TDHamiltonian(time, false);
    auto& Output_ = Output.get_Operator(Space::k);
    auto& Input_ = Input.get_Operator(Space::k);
    auto& H_ = H.get_Operator(Space::k);
    
    commutator(Output_, -im, H_, Input_, false);
    //assert(Output_.is_hermitian());
};