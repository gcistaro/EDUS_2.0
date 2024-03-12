#include <math.h>
#include "Constants.hpp"
#include "Geometry/Coordinate.hpp"
#include "Operator/Operator.hpp"
#include "fftPair/fftPair.hpp"

//computes H_v(x)
double struve(const double& x, const double& v)
{
    if(abs(x)<threshold){
        return 0.;
    }

    double H = 0;
    for(int k=0; k<200; k++){
        int sign = (k%2 == 0 ? 1 : -1);
        H += sign*pow(.5*x, 2*k+v+1)/(std::tgamma(k+1.5)*std::tgamma(k+v+1.5));
    }
    return H;
}

//computes Y_v(x)


class Coulomb
{
    private:
        enum Dimensionality {twoD, threeD} dim;
        double r0=10.;
        double epsilon = 2.;
        double Area=1.;

        Matrix<std::complex<double>> RytovaKeldish_TB;//it contains W[n]=W(rn)
    public:
        Coulomb(const std::array<Operator<std::complex<double>>,3>& r, const int& dim_);
        std::complex<double> W(const Coordinate<k>& q);
        std::complex<double> W(const Coordinate<R>& r);


};


Coulomb::Coulomb(const std::array<Operator<std::complex<double>>,3>& r, const int& dim_)//Pass rR[R=0]
{
    if(dim_ != 2){
        std::cout << "Warning ! Only 2d coulomb is implemented!\n";
        exit(1);
    }
    
    if(dim_==2) dim=twoD;
    if(dim_==3) dim=threeD;

    auto& x = (r[0].get_Operator_R())[Coordinate<R>(0,0,0)];
    auto& y = (r[1].get_Operator_R())[Coordinate<R>(0,0,0)];
    auto& z = (r[2].get_Operator_R())[Coordinate<R>(0,0,0)];
    
    //I use W(n'R1', m'R2', nR1, mR2) = 1/Omega^2*delta(nn')delta(R1R1')delta(mm')delta(R2R2')W(rn+R1-rm-R2)
    RytovaKeldish_TB.initialize(x.get_nrows(), x.get_ncols());
    for(int in=0; in<RytovaKeldish_TB.get_nrows(); in++){
        for(int im=0; im<RytovaKeldish_TB.get_ncols(); im++){
            auto n_atom = Coordinate<R>((r[0].get_Operator_R())[Coordinate<R>(0,0,0)](in,in).real(),
                                               (r[1].get_Operator_R())[Coordinate<R>(0,0,0)](in,in).real(),
                                               (r[2].get_Operator_R())[Coordinate<R>(0,0,0)](in,in).real());
            auto m_atom = Coordinate<R>((r[0].get_Operator_R())[Coordinate<R>(0,0,0)](in,in).real(),
                                               (r[1].get_Operator_R())[Coordinate<R>(0,0,0)](in,in).real(),
                                               (r[2].get_Operator_R())[Coordinate<R>(0,0,0)](in,in).real());

            RytovaKeldish_TB(in, im) = W(n_atom - m_atom);
        }
    }
}

std::complex<double> Coulomb::W(const Coordinate<k>& q)
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

std::complex<double> Coulomb::W(const Coordinate<R>& r)
{
    std::complex<double> Wr;
    auto r_norm = r.norm();
    switch(dim)
    {
        case twoD:
        {
            if(r_norm < threshold){
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
    }
    return Wr;
}



//Hartree: H_{nn'}(R1) = \sum_m \delta{nn'} W(rn-rm+R1) \rho_{mm}(0)

//Fock: H_{nn'}(R1) = 

