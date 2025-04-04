#include "ModelCoulomb.hpp"

/// @brief Calculates the struve special function using the algorithm from XATU
/// @param X__ variable on which we calculate the struve function
/// @param v__ Order of the struve function; only 0 is allowed
/// @return H_v(X)
double struve(const double& X__, const double& v__)
{
    if (v__!=0) {
        throw std::runtime_error("struve functions are defined only for oder 0.");
    }
    double SH0;
    double A0,BY0,P0,Q0,R,S,T,T2,TA0;
	int K, KM;

    S=1.0;
    R=1.0;
    if (X__ <= 20.0) {
        A0=2.0*X__/pi;
        for (K=1; K<61; K++) {
            R=-R*X__/(2.0*K+1.0)*X__/(2.0*K+1.0);
            S=S+R;
            if (fabs(R) < fabs(S)*1.0e-12) goto e15;
        }
        e15:       SH0=A0*S;
    }
    else {
        KM=int(0.5*(X__+1.0));
        if (X__ >= 50.0) KM=25;
        for (K=1; K<=KM; K++) {
           R=-R*pow((2.0*K-1.0)/X__,2);
           S=S+R;
           if (fabs(R) < fabs(S)*1.0e-12) goto e25;
        }
        e25:       T=4.0/X__;
        T2=T*T;
        P0=((((-.37043e-5*T2+.173565e-4)*T2-.487613e-4)*T2+.17343e-3)*T2-0.1753062e-2)*T2+.3989422793;
        Q0=T*(((((.32312e-5*T2-0.142078e-4)*T2+0.342468e-4)*T2-0.869791e-4)*T2+0.4564324e-3)*T2-0.0124669441);
        TA0=X__-0.25*pi;
        BY0=2.0/sqrt(X__)*(P0*sin(TA0)+Q0*cos(TA0));
        SH0=2.0/(pi*X__)*S+BY0;
    }    
    return SH0;
}

/// @brief Constructor that calls the initialization of the class 
/// @param r Position operator read from Wannier90, in a.u.
/// @param dim_ Dimension of the system: 2->monolayer; 3->bulk
/// @param MasterRGrid Grid in R space used in the code, (N.B: it is different than the one in r because that one is read from _tb file)
ModelCoulomb::ModelCoulomb(const std::array<Operator<std::complex<double>>,3>& r, const int& dim_,            
                             const std::shared_ptr<MeshGrid>& MasterRGrid)
{
    initialize(r, dim_, MasterRGrid);
}

/// @brief Initialization of both bare and screened matrix elements calling back the right function
/// @param Rgrid__ MeshGrid where the potential is defined. Note: it must be gamma centered
/// @param nbnd__ Number of bands of the simulation
/// @param Potential__ Matrix elements of the potential
/// @param wannier_centers__ Centers of the wannier functions
/// @param bare__ Boolean: true -> initialize bare, false -> initialize screened
void ModelCoulomb::initialize_Potential( const std::shared_ptr<MeshGrid>& Rgrid__, const int& nbnd__, mdarray<std::complex<double>,3>& Potential__, 
                               const std::vector<Coordinate>& wannier_centers__, const bool& bare__ )
{
    auto size_MG_global =  Rgrid__->get_TotalSize();
    auto size_MG_local = Rgrid__->get_LocalSize();
    Potential__.initialize({size_MG_local, nbnd__, nbnd__});

    #pragma omp parallel for
    for( int iR_local = 0; iR_local < size_MG_local; ++iR_local ) {
        
        auto iR_global = int( Rgrid__->mpindex.loc1D_to_glob1D(iR_local) );
        auto& R = (*Rgrid__)[iR_global];

        for(int in=0; in<ScreenedPotential_.get_Size(1); in++) {
            for(int im=0; im<ScreenedPotential_.get_Size(2); im++){
                auto r = wannier_centers__[in] - R - wannier_centers__[im];
                Potential__(iR_local, in, im) = (bare__ ? V(r) : W(r));
            }
        }
    }
}

/// @brief Initialization of the variables of the class
/// @param r__ The position operator read from Wannier90, in a.u.
/// @param dim__ The dimension of the system: 2->monolayer; 3->bulk
/// @param MasterRGrid__ The grid in R space used in the code, (N.B: it is different than the one in r because that one is read from _tb file)
void ModelCoulomb::initialize(const std::array<Operator<std::complex<double>>,3>& r__, const int& dim__,            
                             const std::shared_ptr<MeshGrid>& MasterRGrid__)
{
    
    if(dim__ != 2){
        throw std::runtime_error("Warning ! Only 2d coulomb is implemented!");
    }
    
    if( dim__ == 2 ) dim_ = twoD;
    if( dim__ == 3 ) dim_ = threeD;

    /* define wannier centers */
    auto index_origin = r__[0].get_Operator_R().get_MeshGrid()->find(Coordinate(0,0,0));
    auto& x0 = (r__[0].get_Operator_R())[index_origin];
    auto& y0 = (r__[1].get_Operator_R())[index_origin];
    auto& z0 = (r__[2].get_Operator_R())[index_origin];

    std::vector<Coordinate> wannier_centers(ScreenedPotential_.get_Size(1));
    for(int in=0; in<ScreenedPotential_.get_Size(1); in++) {
        wannier_centers.push_back( Coordinate(x0(in,in).real(), y0(in,in).real(), z0(in,in).real()));
    }

    /* Define grid centered at 0 */
    Rgrid_ = std::make_shared<MeshGrid>(get_GammaCentered_grid(*MasterRGrid__));

    /* initialize screened and bare potentials matrix elements */
    initialize_Potential( Rgrid_, r__[0].get_Operator(R).get_nrows(), BarePotential_    , wannier_centers, true );
    initialize_Potential( Rgrid_, r__[0].get_Operator(R).get_nrows(), ScreenedPotential_, wannier_centers, false);
}

/// @brief Calculation of the screened interaction on a point r in real space
/// @param r The point on which we calculate the screened interaction
/// @return The screened interaction on the point r: @f$ W(r) = \frac{\pi e^2}{4\pi \epsilon_0 \epsilon_r r_0} \Big[ H_0\Big(\frac{|r|}{r_0}\Big) -Y_0\Big(\frac{|r|}{r_0}\Big)\Big]
std::complex<double> ModelCoulomb::W(const Coordinate& r__)
{
    std::complex<double> Wr;
    auto r_norm = r__.norm();
    switch(dim_)
    {
        case twoD:
        {
            Wr = ( (r_norm < threshold) ? 0. : pi/(r0_*epsilon_)*(struve(r_norm/r0_,0)-y0(r_norm/r0_)));
            break;
        }
        case threeD:
        {
            throw std::runtime_error("Not implemented!");
            break;
        }
        default:
            break;
    }
    return Wr;
}

/// @brief Calculation of the bare interaction on a point r in real space
/// @param r The point on which we calculate the bare interaction
/// @return The bare interaction on the point r: @f$ V(r) = 2./|r| @f$. The factor 2 is for spin degeneracy
std::complex<double> ModelCoulomb::V(const Coordinate& r)
{
    std::complex<double> Vr;
    auto r_norm = r.norm();
    
    Vr = ( (r_norm < threshold) ? 0. : 2./r_norm );
    
    return Vr;
}

/// @brief Set Epsilon (macroscopic dielectric constant)
/// @param Epsilon__ Value we want to use as epsilon
void ModelCoulomb::set_epsilon(const bool& Epsilon__)
{
    epsilon_ =Epsilon__;
}

/// @brief Set r0 in RytovaKeldysh model, not used otherwise
/// @param r0__ Value we want to use as r0 (in a.u.)
void ModelCoulomb::set_r0(const bool& r0__)
{
    r0_ = r0__;
}
