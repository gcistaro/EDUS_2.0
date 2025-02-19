#include "mdContainers/mdContainers.hpp"
#include "RungeKutta/RungeKutta.hpp"
#include "Laser/Laser.hpp"
/* this program can propagate the Newton EOM 
        F=mdv/dt
        v=dr/dt

        with initial condition 
        r0=0
        v0=0
 */

int NumberOfParticles = 1;
double InitialTime = 0.;
double koverm = 0.01;
double ResolutionTime = 0.1;





int main()
{
    mdarray<double, 3> XandV( {2, NumberOfParticles, 3} ); //first index  -> {x, \dot{x}}
                                                      //second index -> Particle
                                                      //third index  -> cartesian direction 

    mdarray<double,2> X(&XandV(0,0,0), {NumberOfParticles, 3});
    mdarray<double,2> V(&XandV(1,0,0), {NumberOfParticles, 3});
    
    
    Laser laser;
    laser.set_Intensity(1.e+06, Wcm2);
    laser.set_Omega(2*3.14*std::sqrt(0.01), AuEnergy);
    laser.set_Polarization(Coordinate(1,0,0));
    laser.set_NumberOfCycles(1);

    auto Force = [&](const auto& X, const auto& V, const double& time){
    mdarray<double, 2> force( { NumberOfParticles, 3 });
    force.fill(0.);
    /* internal forces */
    for( int iparticle = 0; iparticle < NumberOfParticles; ++iparticle ) {
        for( int ix : {0, 1, 2} ) {
            force( iparticle, ix ) += -koverm*X( iparticle, ix );
        }
    }

    /* external forces */
    for( int iparticle = 0; iparticle < NumberOfParticles; ++iparticle ) {
        for( int ix : {0, 1, 2} ) {
            force( iparticle, ix ) += laser(time).get("Cartesian")[0];
        }
    }

    return force;
};



    std::function<void(mdarray<double, 3>&)>  
    InitialCondition = [&](auto& Function){
        Function.fill(0.);
        Function(0,0,0) = 0.;
    };

    std::function<void(mdarray<double, 3>&, double const&, mdarray<double, 3> const&)>
    SourceTerm = [&](auto& Output, const double InputTime, const auto& InputFunction){
        mdarray<double,2> X0(&(const_cast<mdarray<double, 3>&>(InputFunction))(0,0,0), {NumberOfParticles, 3});
        mdarray<double,2> V0(&(const_cast<mdarray<double, 3>&>(InputFunction))(1,0,0), {NumberOfParticles, 3});

        
        mdarray<double,2> X1(&(const_cast<mdarray<double, 3>&>(Output))(0,0,0), {NumberOfParticles, 3});
        mdarray<double,2> V1(&(const_cast<mdarray<double, 3>&>(Output))(1,0,0), {NumberOfParticles, 3});


        /* dr/dt = v */
        std::copy( V0.begin(), V0.end(), X1.begin() );

        /* dv/dt = F */
        auto force = Force(X0, V0, InputTime);
        std::copy(force.begin(), force.end(), V1.begin());
    };


    auto rungekutta = RungeKutta<mdarray<double,3>>(XandV, InitialCondition, SourceTerm);
    rungekutta.set_InitialTime(InitialTime);
    rungekutta.set_ResolutionTime(ResolutionTime);


    std::ofstream fp_out("t_X_V.txt");
    std::ofstream fp_Force("Force.txt");
    std::ofstream fp_extForce("ExternalForce.txt");

    for(double it=0; it<=10000; it++){
        auto time = rungekutta.get_CurrentTime();
        fp_Force << time << " " << Force(X,V,time)(0,0) << std::endl;
        fp_out << time << " " << X(0,0) << " " << V(0,0) << std::endl;
        fp_extForce << time << " " << laser(time).get("Cartesian")[0] << std::endl;
        rungekutta.Propagate();
    }

    fp_out.close();
    fp_Force.close();
    fp_extForce.close();
}