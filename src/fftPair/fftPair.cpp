#include "fftPair/fftPair.hpp"
#include <algorithm> 
#include <functional> 
#include <cctype>
#include <locale>
#include "core/profiler.hpp"

FourierTransform::FourierTransform
(mdarray<std::complex<double>, 2>& Array_x__, mdarray<std::complex<double>, 2>& Array_k__, const std::vector<int>& Dimensions__)
{
    this->initialize(Array_x__, Array_k__, Dimensions__);
}

void 
FourierTransform::initialize
(mdarray<std::complex<double>, 2>& Array_x__, mdarray<std::complex<double>, 2>& Array_k__, const std::vector<int>& Dimensions__, 
            const std::string tagname_)
{
    tagname = tagname_;
    IsFFT = true;
    Array_x = &Array_x__;
    Array_k = &Array_k__;

#ifdef EDUS_MPI
    howmany = Array_x->get_Size(1);
#else
    howmany = Array_x->get_Size(0);
#endif
    TotalSize = 1;
    for(auto& dim__ : Dimensions__) {
        TotalSize *= dim__; 
    }
    SqrtTotalSize = std::sqrt(double(TotalSize));
    dim = Dimensions__.size();
    Dimensions = Dimensions__;
    idist = TotalSize;
    odist = TotalSize;

//initializing plans, two for each object for +1 and -1 transforms.
#ifdef EDUS_MPI
    std::ptrdiff_t* Dimensions_ptr = new std::ptrdiff_t[Dimensions.size()];
    for( auto idim = 0; idim < int(Dimensions.size()); ++idim ) {
        Dimensions_ptr[idim] = Dimensions[idim];
    }
    //auxiliary variable needed 
    //from x to k (fft to Fourier space)
    MyPlan_FWD = fftw_mpi_plan_many_dft(dim, Dimensions_ptr,
                                 howmany, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK,
                                 reinterpret_cast<fftw_complex*>(&(*Array_x)[0]),
                                 reinterpret_cast<fftw_complex*>(&(*Array_k)[0]), 
                                 MPI_COMM_WORLD, -1, FFTW_ESTIMATE);

    //from k to x (fft to original space) -> sign = +1 correspond to ifft -> f(x) = sum_n c_n e^{+inx}
    MyPlan_BWD = fftw_mpi_plan_many_dft(dim, Dimensions_ptr,
                                 howmany, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK,
                                 reinterpret_cast<fftw_complex*>(&(*Array_k)[0]),
                                 reinterpret_cast<fftw_complex*>(&(*Array_x)[0]), 
                                 MPI_COMM_WORLD, +1, FFTW_ESTIMATE);
    //as per user guide: this corresponds to the same parameter in the serial advanced interface 
    //(see Advanced Complex DFTs) with stride = howmany and dist = 1. Meaning that data are contiguous 
    //in the dimension of the bands (no transpose needed!)
    delete[] Dimensions_ptr;
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
    IsFFT = false;
    Mesh = Mesh__;
    howmany = Array_x__.get_Size(0);
    dim = Mesh[0].size();
    TotalSize = Array_x__.get_Size(1);
    SqrtTotalSize = std::sqrt(double(TotalSize));
    Array_x = &Array_x__;
}

void FourierTransform::fft(const int& sign)
{
    assert(IsFFT);
    assert( sign == 1 || sign == -1 );
    auto& output = (sign == +1 ? (*Array_x) : (*Array_k) ); 
    auto& MyPlan = (sign == +1 ? (MyPlan_BWD) : (MyPlan_FWD) );

    fftw_execute(MyPlan);

    if( sign == -1 ) {
        #pragma omp parallel for schedule(static)
        for(int index = 0; index < output.get_TotalSize(); ++index) {
            auto& output_el = output[index];
            output_el /= TotalSize;
        }
    }
}


std::complex<double> FourierTransform::dft(const std::vector<double>& Point, const int& h, const int& sign) 
{
    assert(int(Point.size()) == dim);
    //mdarray<std::complex<double>, 1> FT({Array_x->get_Size()[0]});
    std::complex<double> FT = 0.;//.fill(std::complex<double>(0.));

    //std::complex<double> FourierTransform = 0;
    static std::complex<double> im2pi = im*2.*pi;
    double DotProduct;

        for(int i=0; i<int(Mesh.size()); i++){
            DotProduct = 0;
            for(int ix=0; ix<dim; ix++){
                DotProduct += Point[ix]*Mesh[i][ix];
            }
            FT += std::exp(double(sign)*im2pi*DotProduct)*(*Array_x)(h,i);
	    }
    return FT;
}

mdarray<std::complex<double>, 2> FourierTransform::dft(const std::vector<std::vector<double>>& ArrayOfPoints, const int& sign)
{
    PROFILE("dft");
    auto md_K = mdarray<std::complex<double>, 2>({Array_x->get_Size()[0], int(ArrayOfPoints.size())});
    Array_k = &md_K;
    #pragma omp parallel for collapse(2)
    for(int h=0; h<howmany; ++h){
        for(int ip=0; ip<int(ArrayOfPoints.size()); ++ip){
            (*Array_k)(h,ip) = dft(ArrayOfPoints[ip], h, sign);
        }
    }
    return (*Array_k);
}


FourierTransform::~FourierTransform()
{
    if( IsFFT && destruct )
     {
        fftw_destroy_plan(MyPlan_FWD);
        fftw_destroy_plan(MyPlan_BWD);
        destruct = false;
    }
}
