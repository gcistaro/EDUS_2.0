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
    /* next line aligns k and R components of input (DM_{i-1}) */
    const_cast<Operator<std::complex<double>>&>(Input__).go_to_R(true);
    
    /* Gradient term:    Output +=   (E.Nabla) * Input */ 
    Output__.lock_space(SpaceOfPropagation_Gradient_);
    if (!ctx_->cfg().peierls()) {
        kgradient_.Calculate(1.+0.*im, Output__.get_Operator(SpaceOfPropagation_Gradient_), 
                        Input__.get_Operator(SpaceOfPropagation_Gradient_), 
                        setoflaser_(time__), true);

    }

    /* IPA Hamiltonian H_ = H0_ + E \cdot r*/
    Calculate_TDHamiltonian(time__, true);
    
    /* Peierls transformation H_(R) = H_(R)*exp(+i*A(t) \cdot R) */
    if(ctx_->cfg().peierls()) {
        apply_peierls_phase(H_, time__, +1);
    }

    /* Coulomb interaction H_ += \Sigma^H + \Sigma^{SEX} */
    H_.go_to_R();
    coulomb_.EffectiveHamiltonian( H_, Input__, false); 


    // Output__ += -i * [ H_, Input__ ]
    Output__.go_to_k();
    H_.go_to_k();
    const_cast<Operator<std::complex<double>>&>(Input__).lock_space(Space::k); //this is already updated, no need to ft

    auto& Output = Output__.get_Operator(SpaceOfPropagation_);
    auto& Input = Input__.get_Operator(SpaceOfPropagation_);
    auto& H = H_.get_Operator(SpaceOfPropagation_);
    commutator(Output, -im, H, Input, false);
};