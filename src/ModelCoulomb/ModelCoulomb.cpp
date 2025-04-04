#include "ModelCoulomb.hpp"

/// @brief Calculates the struve special function using the algorithm from XATU
/// @param X variable on which we calculate the struve function
/// @param v Order of the struve function; only 0 is allowed
/// @return H_v(X)
double struve(const double& X, const double& v)
{
    if (v!=0) {
        throw std::runtime_error("struve functions are defined only for oder 0.");
    }
    double SH0;
    double A0,BY0,P0,Q0,R,S,T,T2,TA0;
	int K, KM;

    S=1.0;
    R=1.0;
    if (X <= 20.0) {
        A0=2.0*X/pi;
        for (K=1; K<61; K++) {
            R=-R*X/(2.0*K+1.0)*X/(2.0*K+1.0);
            S=S+R;
            if (fabs(R) < fabs(S)*1.0e-12) goto e15;
        }
        e15:       SH0=A0*S;
    }
    else {
        KM=int(0.5*(X+1.0));
        if (X >= 50.0) KM=25;
        for (K=1; K<=KM; K++) {
           R=-R*pow((2.0*K-1.0)/X,2);
           S=S+R;
           if (fabs(R) < fabs(S)*1.0e-12) goto e25;
        }
        e25:       T=4.0/X;
        T2=T*T;
        P0=((((-.37043e-5*T2+.173565e-4)*T2-.487613e-4)*T2+.17343e-3)*T2-0.1753062e-2)*T2+.3989422793;
        Q0=T*(((((.32312e-5*T2-0.142078e-4)*T2+0.342468e-4)*T2-0.869791e-4)*T2+0.4564324e-3)*T2-0.0124669441);
        TA0=X-0.25*pi;
        BY0=2.0/sqrt(X)*(P0*sin(TA0)+Q0*cos(TA0));
        SH0=2.0/(pi*X)*S+BY0;
    }    
    return SH0;
}

/// @brief Constructor that calls the initialization of the class 
/// @param r The position operator read from Wannier90, in a.u.
/// @param dim_ The dimension of the system: 2->monolayer; 3->bulk
/// @param MasterRGrid The grid in R space used in the code, (N.B: it is different than the one in r because that one is read from _tb file)
ModelCoulomb::ModelCoulomb(const std::array<Operator<std::complex<double>>,3>& r, const int& dim_,            
                             const std::shared_ptr<MeshGrid>& MasterRGrid)
{
    initialize(r, dim_, MasterRGrid);
}

/// @brief Initialization of the variables of the class
/// @param r The position operator read from Wannier90, in a.u.
/// @param dim_ The dimension of the system: 2->monolayer; 3->bulk
/// @param MasterRGrid The grid in R space used in the code, (N.B: it is different than the one in r because that one is read from _tb file)
void ModelCoulomb::initialize(const std::array<Operator<std::complex<double>>,3>& r, const int& dim_,            
                             const std::shared_ptr<MeshGrid>& MasterRGrid)
{
    
    if(dim_ != 2){
        throw std::runtime_error("Warning ! Only 2d coulomb is implemented!");
    }
    
    if(dim_==2) dim=twoD;
    if(dim_==3) dim=threeD;
    auto index_origin = r[0].get_Operator_R().get_MeshGrid()->find(Coordinate(0,0,0));

    auto& x0 = (r[0].get_Operator_R())[index_origin];
    auto& y0 = (r[1].get_Operator_R())[index_origin];
    auto& z0 = (r[2].get_Operator_R())[index_origin];

    //I use W(n'R1', m'R2', nR1, mR2) = 1/Omega^2*delta(nn')delta(R1R1')delta(mm')delta(R2R2')W(rn+R1-rm-R2)
    ScreenedPotential_.initialize({MasterRGrid->get_TotalSize(),
                                 r[0].get_Operator_R().get_nrows(), 
                                 r[0].get_Operator_R().get_ncols()});
    Rgrid = std::make_shared<MeshGrid>(get_GammaCentered_grid(*MasterRGrid));

//    #pragma omp parallel for schedule(static)
    for(int iR=0; iR<Rgrid->get_TotalSize(); ++iR) {//TB.get_Size(0); ++iR) {
        auto& Rcart = (*Rgrid)[iR].get("Cartesian");
        for(int in=0; in<ScreenedPotential_.get_Size(1); in++){
            //ratom_n = rn - R
            auto ratom_n = Coordinate(x0(in,in).real() - Rcart[0],
                                      y0(in,in).real() - Rcart[1],
                                      z0(in,in).real() - Rcart[2]);

            for(int im=0; im<ScreenedPotential_.get_Size(2); im++){
                //ratom_m = rm + S
                auto ratom_m = Coordinate(x0(im,im).real(),//+Scart[0],
                                          y0(im,im).real(),//+Scart[1],
                                          z0(im,im).real());//+Scart[2]);
                ScreenedPotential_(iR, in, im) = Potential(ratom_n - ratom_m);
            }
        }
    }
}

/// @brief Calculation of the screened interaction on a point r in real space
/// @param r The point on which we calculate the screened interaction
/// @return The screened interaction on the point r: W(r)
std::complex<double> ModelCoulomb::W(const Coordinate& r)
{
    std::complex<double> Wr;
    auto r_norm = r.norm();
    switch(dim)
    {
        case twoD:
        {
            if(r_norm < threshold ){
                Wr = 0.;
            }
            else{
                Wr = pi/(r0*epsilon)*(struve(r_norm/r0,0)-y0(r_norm/r0));
            }
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
/// @return The bare interaction on the point r: @f$ V(r) = 2./r @f$. The factor 2 is for spin degeneracy
std::complex<double> ModelCoulomb::V(const Coordinate& r)
{
    std::complex<double> Vr;
    auto r_norm = r.norm();
    
    if(r_norm < threshold ){
        Vr = 0.;
    }
    else{
        Vr = 2./r_norm;
    }

    return Vr;
}

