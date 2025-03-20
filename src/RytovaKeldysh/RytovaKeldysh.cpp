#include "RytovaKeldysh.hpp"

//computes H_v(x)
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

RytovaKeldysh::RytovaKeldysh(const std::array<Operator<std::complex<double>>,3>& r, const int& dim_,            
                             const std::shared_ptr<MeshGrid>& MasterRGrid)
{
    initialize(r, dim_, MasterRGrid);
}

void RytovaKeldysh::initialize(const std::array<Operator<std::complex<double>>,3>& r, const int& dim_,            
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
    TB.initialize({MasterRGrid->get_TotalSize(),
                                 r[0].get_Operator_R().get_nrows(), 
                                 r[0].get_Operator_R().get_ncols()});
    Rgrid = std::make_shared<MeshGrid>(get_GammaCentered_grid(*MasterRGrid));

//    #pragma omp parallel for schedule(static)
    for(int iR=0; iR<Rgrid->get_TotalSize(); ++iR) {//TB.get_Size(0); ++iR) {
        auto& Rcart = (*Rgrid)[iR].get("Cartesian");
        for(int in=0; in<TB.get_Size(1); in++){
            //ratom_n = rn - R
            auto ratom_n = Coordinate(x0(in,in).real() - Rcart[0],
                                      y0(in,in).real() - Rcart[1],
                                      z0(in,in).real() - Rcart[2]);

            for(int im=0; im<TB.get_Size(2); im++){
                //ratom_m = rm + S
                auto ratom_m = Coordinate(x0(im,im).real(),//+Scart[0],
                                          y0(im,im).real(),//+Scart[1],
                                          z0(im,im).real());//+Scart[2]);
                TB(iR, in, im) = Potential(ratom_n - ratom_m);
            }
        }
    }
}

/*
std::complex<double> Coulomb::W(const Coordinate& q)
{

    double qnorm = q.norm();
    std::complex<double> Wq;
    switch(dim)
    {
        case twoD:
        {
            //Rytova-Keldish potential
            Wq = 1./(epsilon*Area*qnorm*(r0*qnorm+1.));
            break;
        }
        case threeD:
        {
            throw std::runtime_error("threed not implemented");
            break;
        }
    }
    return Wq;
}
*/
std::complex<double> RytovaKeldysh::Potential(const Coordinate& r)
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

            //done via Fourier transform dft on W(q)
            //set up W(q)
/*            MeshGrid<k> aux_mg(50., {50, 50, 50});
            mdarray<std::complex<double>, 2> Wq({1, aux_mg.get_mesh().size()});
            for(int im=0; im<aux_mg.get_mesh().size(); im++){
                Wq(0,im) = W(aux_mg[im]);
            }

            //go to plain MeshGrid and plain r
            std::vector<std::vector<double>> plain_mg(aux_mg.get_mesh().size(), std::vector<double>(3));
            for(int im=0; im<plain_mg.size(); im++){
                auto& aux_vec = aux_mg[im].get("LatticeVectors");
                for(auto& ix : {0,1,2}){
                    plain_mg[im][ix] = aux_vec[ix];
                }
            }
            auto& aux_r = r.get("LatticeVectors");
            std::vector<double> plain_r(3);
            for(auto& ix : {0,1,2}){
                plain_r[ix] = aux_r(ix);
            }

            //dft on W(q)
            FourierTransform ft_;
            ft_.initialize(Wq, plain_mg);

            Wr = ft_.dft(plain_r, +1)(0);
*/
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

