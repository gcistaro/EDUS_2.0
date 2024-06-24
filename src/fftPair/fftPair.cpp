#include "fftPair/fftPair.hpp"

FourierTransform::FourierTransform
(mdarray<std::complex<double>, 2>& Array_x__, mdarray<std::complex<double>, 2>& Array_k__, const std::vector<int>& Dimensions__)
{
    this->initialize(Array_x__, Array_k__, Dimensions__);
}

void 
FourierTransform::initialize
(mdarray<std::complex<double>, 2>& Array_x__, mdarray<std::complex<double>, 2>& Array_k__, const std::vector<int>& Dimensions__)
{
    Array_x = &Array_x__;
    Array_k = &Array_k__;
    howmany = Array_x->get_Size(0);
    TotalSize = Array_x->get_Size(1);
    SqrtTotalSize = std::sqrt(double(TotalSize));
    dim = Dimensions__.size();
    Dimensions = Dimensions__;
    idist = TotalSize;
    odist = TotalSize;

//initializing plans, two for each object for +1 and -1 transforms.
#ifdef NEGF_MPI
//TODO!! HERE YOU NEED TO ADD THE CONSTRUCTION OF THE PLAN!!!!!
#else
    //from x to k (fft to Fourier space)
    MyPlan_FWD = fftw_plan_many_dft(dim, &Dimensions[0], howmany,
                                    reinterpret_cast<fftw_complex*>(&(*Array_x)[0]), inembed, istride, idist, 
                                    reinterpret_cast<fftw_complex*>(&(*Array_k)[0]), onembed, ostride, odist,
                                    -1, FFTW_ESTIMATE);
    //from k to x (fft to original space) -> sign = +1 correspond to ifft -> f(x) = sum_n c_n e^{+inx}
    MyPlan_BWD = fftw_plan_many_dft(dim, &Dimensions[0], howmany,
                                    reinterpret_cast<fftw_complex*>(&(*Array_k)[0]), inembed, istride, idist, 
                                    reinterpret_cast<fftw_complex*>(&(*Array_x)[0]), onembed, ostride, odist,
                                    +1, FFTW_ESTIMATE);
#endif
}

void 
FourierTransform::initialize
(mdarray<std::complex<double>, 2>& Array_x__, const std::vector<std::vector<double>>& Mesh__)
{
    Mesh = Mesh__;
    howmany = Array_x__.get_Size(0);
    dim = Mesh[0].size();
    TotalSize = Array_x__.get_Size(1);
    SqrtTotalSize = std::sqrt(double(TotalSize));
    Array_x = &Array_x__;
}

void FourierTransform::fft(const int& sign)
{
    auto& output = (sign == +1 ? (*Array_x) : (*Array_k) ); 
    auto& MyPlan = (sign == +1 ? (MyPlan_BWD) : (MyPlan_FWD) );

    fftw_execute(MyPlan);

    #pragma omp parallel for schedule(static)
    for(auto& output_el : output){
        output_el /= SqrtTotalSize;
    }
}


mdarray<std::complex<double>, 1> FourierTransform::dft(const std::vector<double>& Point, const int& sign) 
{
    assert(Point.size() == dim);
    mdarray<std::complex<double>, 1> FT({Array_x->get_Size()[0]});
    FT.fill(std::complex<double>(0.));

    //std::complex<double> FourierTransform = 0;
    static std::complex<double> im2pi = im*2.*pi;
    for(int h=0; h<howmany; h++){
        for(int i=0; i<Mesh.size(); i++){
            double DotProduct = 0;
            for(int ix=0; ix<dim; ix++){
                DotProduct += Point[ix]*Mesh[i][ix];
            }
            FT(h) += std::exp(double(sign)*im2pi*DotProduct)*(*Array_x)(h,i);
	    }
    }
    return FT;
}

mdarray<std::complex<double>, 2> FourierTransform::dft(const std::vector<std::vector<double>>& ArrayOfPoints, const int& sign)
{
    auto md_K = mdarray<std::complex<double>, 2>({Array_x->get_Size()[0], ArrayOfPoints.size()});
    Array_k = &md_K;
    for(int ip=0; ip<ArrayOfPoints.size(); ++ip){
        auto FT = dft(ArrayOfPoints[ip], sign);
        //copy to Array_k
        for(int h=0; h<howmany; ++h){
            (*Array_k)(h,ip) = FT(h);
	    }
    }
    return (*Array_k);
}


FourierTransform::~FourierTransform()
{
    fftw_destroy_plan(MyPlan_FWD);
    fftw_destroy_plan(MyPlan_BWD);
}