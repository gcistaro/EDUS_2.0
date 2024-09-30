std::function<void(Operator<std::complex<double>>&, double const&, Operator<std::complex<double>> const&)> 
SourceTerm = 
[&](Operator<std::complex<double>>& Output, const double& time, const Operator<std::complex<double>>& Input)
{
    // Output = -i * [ H, Input ]
    const_cast<Operator<std::complex<double>>&>(Input).get_Operator_k().make_hermitian();
    Calculate_TDHamiltonian(time, true);
    Output.go_to_k(false);
    const_cast<Operator<std::complex<double>>&>(Input).go_to_k(false);
    H.go_to_k(false);
    auto& Output_ = Output.get_Operator(Space::k);
    auto& Input_ = Input.get_Operator(Space::k);
    auto& H_ = H.get_Operator(Space::k);
    
    commutator(Output_, -im, H_, Input_, true);
    assert(Output_.is_hermitian());
    Output.go_to_R();
    const_cast<Operator<std::complex<double>>&>(Input).go_to_R();
    
    //Output +=   (E.Nabla) * Input
    PROFILE_START("i(E.R)*Input");
    kgradient.Calculate(1.+0.*im,Output.get_Operator(Space::R), 
                        Input.get_Operator(Space::R), 
                        laser(time), false);
    PROFILE_STOP("i(E.R)*Input");
    const_cast<Operator<std::complex<double>>&>(Input).go_to_k(false);
    Output.go_to_k();
};