#include "GreenFunction.hpp"

void GreenFunction::initialize( DESolver<Operator<std::complex<double>>>& DEsolver_DM__ )
{
#ifndef EDUS_HDF5
    std::runtime_error("Error in GreenFunction::initialize: Compiled without HDF5 support");
#endif
    /* allocate memory for operators */
    DEsolver_DM_ = &DEsolver_DM__;
    DM_           .initialize_fft(DEsolver_DM__.get_Function());
    H_            .initialize_fft(DEsolver_DM__.get_Function());
    Ut_           .initialize_fft(DEsolver_DM__.get_Function());
    Utime_        .initialize_fft(DEsolver_DM__.get_Function());
    Uptime_       .initialize_fft(DEsolver_DM__.get_Function());
    GR_           .initialize_fft(DEsolver_DM__.get_Function());
    
    /* initialize solver for Ut */
    std::function<void(Operator<std::complex<double>>&)> 
    InitialCondition = 
    [&](Operator<std::complex<double>>& U)
    {
        /* the initial evolution operator is the identity */
        U.get_Operator(Space::k).identity();
    };

    std::function<void(Operator<std::complex<double>>&, double const&, Operator<std::complex<double>> const&)> 
    SourceTerm = 
    [&](Operator<std::complex<double>>& Output, const double& time, const Operator<std::complex<double>>& Input)
    {
        /* EOM:      i\dot{U} = H.U */
        /* read H from saved datas */
        auto CurrentTime = DEsolver_DM_->get_CurrentTime();
        std::string name_ = "output.h5";
        auto node = get_it_sparse(CurrentTime);
        H_.get_Operator(Space::k).load(name_, node, nodename::H0pCoulomb); 
        /* Output = -i H.U */
        multiply(Output.get_Operator(Space::k), -im, 
                 H_.get_Operator(Space::k), 
                 Input.get_Operator(Space::k));
    };
    DEsolver_Ut_.initialize( Ut_, InitialCondition, SourceTerm, AB, 5 );
}

void GreenFunction::Propagator(Operator<std::complex<double>>& Ut__, const double t__)
{
    if ( ( DEsolver_Ut_.get_CurrentTime() - t__ ) > 1.e-06 ) {
        //read it from file
    }
    else {
        auto nsteps = (t__ - DEsolver_Ut_.get_CurrentTime())/double(DEsolver_Ut_.get_ResolutionTime() );
        DEsolver_Ut_.Propagate( nsteps );
    }
}

/*
    We calculate 
    G<(t,t') = Theta(t-t')U(t)U^\dagger(t')G<(t',t') + Theta(t'-t) G<(t,t)U(t)U^\dagger(t')
    as formula 31 of Stefanucci, An essential Introduction to NEGF methods for Real-Time simulations
    Note that G<(t',t) = Theta(t'-t)U(t')U^\dagger(t)G<(t,t) + Theta(t-t') G<(t',t')U(t')U^\dagger(t)
                       = Theta(t'-t) ( G<^\dagger(t,t)U(t)U^\dagger(t') )^\dagger 
                         + Theta(t-t') ( U(t)U^\dagger(t')G<^\dagger(t',t') )\dagger
                       = G<(t,t')^\dagger

    We also do the change of variables: 
              trel   = t'-t      |_______     t' = tave + trel/2
                                 |_______>
              tave = (t+t')/2    |            t  = tave - trel/2

*/

BlockMatrix<std::complex<double>>& GreenFunction::GKBA_rel_ave(const double trel__, const double tave__)
{
    return GKBA(tave__+trel__/2., tave__-trel__/2.);
}

BlockMatrix<std::complex<double>>& GreenFunction::GKBA(const double time__, const double ptime__)
{
    auto Theta = [](const double& t, const double& pt) { return ( t - pt ) > 1.e-06;};
    if( !Theta(time__,ptime__) ) {
        BlockMatrix<std::complex<double>>* aux = &GKBA( ptime__,time__ );
        aux->make_dagger();
        return *aux;
    }

    /* DM_ <- rho(t') */
    std::string name_ = "output.h5";
    auto node = get_it_sparse(  ptime__ );
    DM_.get_Operator(Space::k).load(name_, node, nodename::DMk); 

    if( std::abs(time__-ptime__) < 1.e-06 ) {
        return DM_.get_Operator(Space::k);
    }

    /* here we calculate U(t)U^\dagger(t')G<(t',t') */

    /* Utime_ <- U(t) */
    Propagator(Utime_, time__);
    /* Uptime_ <- U(t') */
    Propagator(Uptime_, ptime__);

    /* GR_ = U(t)U^\dagger(t') */
    multiply( GR_.get_Operator(Space::k), 1.+0.*im, 
              Utime_.get_Operator(Space::k), Uptime_.get_Operator(Space::k));
    /* G< = GR_*DM */
    /* we use Utime_ array as it is not needed anymore */
    auto& Lesser = Utime_.get_Operator(Space::k);
    multiply( Lesser, 1.+0.*im, 
              GR_.get_Operator(Space::k), DM_.get_Operator(Space::k) );
    return Lesser;
}


int GreenFunction::get_it_sparse(const double& time) const
{
    return int(round(time/DEsolver_DM_->get_ResolutionTime()/PrintResolution_));
}
