#include "RytovaKeldysh.hpp"

//computes H_v(x)
double struve(const double& x, const double& v)
{
    if(abs(x)<threshold){
        return 0.;
    }

    double H = 0;
    for(int k=0; k<100; k++){
        auto tgamma__ = std::tgamma(k+v+1.5);
        int sign = (k%2 == 0 ? 1 : -1);
        H += sign*pow(.5*x, 2*k+v+1)/(std::tgamma(k+1.5)*std::tgamma(k+v+1.5));
        //std::cout << "rnorm " << x << " " << sign*pow(.5*x, 2*k+v+1)<<  " k " << k << " v " << v << " H "<< H << " tgamma " << std::tgamma(k+1.5) << " " << std::tgamma(k+v+1.5) << std::endl;
    }
    return H;
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
        std::cout << "Warning ! Only 2d coulomb is implemented!\n";
        exit(1);
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
    Rgrid = MasterRGrid;

    //Calculate values of the potential
    std::cout << "TB.get_Size(0): " << TB.get_Size(0) << std::endl;
    std::cout << "TB.get_Size(1): " << TB.get_Size(1) << std::endl;

//    #pragma omp parallel for schedule(static)
    for(int iR=0; iR<1; ++iR) {//TB.get_Size(0); ++iR) {
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
                std::cout << TB(iR, in, im) << std::endl;
            }
        }
    }
    std::cout << "out of loop!\n";
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
            std::cout << "threed not implemented\n";
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
                Wr = pi/(2.*r0*epsilon)*(struve(r_norm/r0,0)-y0(r_norm/r0));
            }

            //done via Fourier transform dft on W(q)
            //set up W(q)
/*            MeshGrid<k> aux_mg(50., {50, 50, 50});
            for(int im=0; im<aux_mg.get_mesh().size(); im++){
                std::cout << aux_mg[im].get("LatticeVectors") << std::endl;
            }
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
            std::cout << "Not implemented!\n";
            exit(1);
        }
        default:
            break;
    }
        std::cout << r_norm << " "<< Wr << std::endl;

    return Wr;
}

