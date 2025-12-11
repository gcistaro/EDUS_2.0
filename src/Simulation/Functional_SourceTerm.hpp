/// @brief This standard function contains what equals the derivative in time of the density matrix. 
/// From the Schrodinger equation, we know:
/// @f[
/// \frac{\partial \rho}{\partial t} = -i [H_0 + H_{\text{eff}} + \boldsymbol{\varepsilon}(t)\cdot \Xi, \rho] + 
/// \boldsymbol{\varepsilon}(t)\nabla_\textbf{k} \rho 
/// @f]
/// @param Output__ We store here @f$ \frac{\partial \rho}{\partial t} @f$
/// @param time__ Current time of the simulation, to calculate the time dependent hamiltonian and the laser
/// @param Input__ Input density matrix, to be used as the density matrix on the RHS of the equation
std::function<void(Operator<std::complex<double>>&, double const&, Operator<std::complex<double>> const&)> 
SourceTerm = 
[&](Operator<std::complex<double>>& Output__, const double& time__, const Operator<std::complex<double>>& Input__)
{
    
    if ( SpaceOfPropagation_Gradient_ == R ) {
        Output__.go_to_R(false);
        const_cast<Operator<std::complex<double>>&>(Input__).go_to_R(true);
    }

    /* Gradient term:    Output +=   (E.Nabla) * Input 
    it acts only if we don't do Peierls transformation */   
    if (!ctx_->cfg().peierls()) {
        kgradient_.Calculate(1.+0.*im, Output__.get_Operator(SpaceOfPropagation_Gradient_), 
                        Input__.get_Operator(SpaceOfPropagation_Gradient_), 
                        setoflaser_(time__), true);

    }

    /* Coulomb interaction */
    H_.go_to_R(false);
    H_.get_Operator_R().fill(0);
    coulomb_.EffectiveHamiltonian( H_, Input__, true); 

    if ( SpaceOfPropagation_Gradient_ == R ) {
        Output__.go_to_k(true);
        const_cast<Operator<std::complex<double>>&>(Input__).go_to_k(false);
    }


    // Output += -i * [ H, Input ]
    if( SpaceOfCalculateTDHamiltonian_ == k ) H_.go_to_k(true);
    Calculate_TDHamiltonian(time__, false);
    if( SpaceOfCalculateTDHamiltonian_ == R ) H_.go_to_k(true);

    auto& Output = Output__.get_Operator(SpaceOfPropagation_);
    auto& Input = Input__.get_Operator(SpaceOfPropagation_);
    auto& H = H_.get_Operator(SpaceOfPropagation_);
    commutator(Output, -im, H, Input, false);

};