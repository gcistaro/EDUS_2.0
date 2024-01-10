#include <vector>
#include <memory>
#include <complex>
#include "../mdContainers/mdContainers.hpp"
#include <fftw3.h>



template<size_t dim>
class fftPair
{
    private:
        std::shared_ptr<mdarray<std::complex<double>, dim>> InputArray;
        std::shared_ptr<mdarray<std::complex<double>, dim>> OutputArray;
        std::vector<std::vector<double>> Mesh;
        std::vector<std::vector<double>> Mesh_fft;
    	int Dimensions[dim];
        int TotalSize;
        int istride = 1;
        int ostride = 1;
        int idist;
        int odist;
        int* inembed = nullptr;
        int* onembed = nullptr;
        int howmany = 1;

        fftw_plan MyPlan;
    public:
        fftPair(){};
        fftPair(mdarray<std::complex<double>,dim>& Data_)
        {
            Mesh.resize(dim);
            Mesh_fft.resize(dim);
            InputArray = std::make_shared<mdarray<std::complex<double>,dim>>(Data_);
            TotalSize = InputArray->get_TotalSize();
            idist = TotalSize;
            odist = TotalSize;

            for(int i=0; i<dim; i++){
                Dimensions[i] = InputArray->get_Size(i);
                
                Mesh[i].resize(Dimensions[i]);
                Mesh_fft[i].resize(Dimensions[i]);
                //by default, we sample always [0,Npoints) in direct space, i.e. [0,1) in fourier space
                for(int ix=0; ix<Dimensions[i]; ix++){
                    Mesh[i][ix] = double(ix);
                    Mesh_fft[i][ix] = double(ix)/double(Dimensions[i]);
                }
            }
        };

	    void fft(const int& sign);
        std::complex<double> fftPair<dim>::dft(const auto& Point) 
        
        const mdarray<std::complex<double>, dim>& get_OutputArray() const
        {
            return (*OutputArray);
        }

};


template<size_t dim>
void fftPair<dim>::fft(const int& sign)
{
    OutputArray = std::make_shared<mdarray<std::complex<double>, dim>>(mdarray<std::complex<double>, dim>(InputArray->get_Size()));

    //TODO: be sure Data is allocated till the end
    MyPlan = fftw_plan_many_dft(dim, //Dimension of the fft
	                       &Dimensions[0],   
                           howmany,       
                           reinterpret_cast<fftw_complex*>(&(*InputArray)[0]),                    //row-major ordered 
                           inembed,
                           istride,
                           idist,
                           reinterpret_cast<fftw_complex*>(&(*OutputArray)[0]),                   //row-major ordered 
                           onembed,
                           ostride,
                           odist,
                           sign,
			               FFTW_ESTIMATE);
    fftw_execute(MyPlan);
    fftw_destroy_plan(MyPlan);
}




template<size_t dim>
std::complex<double> fftPair<dim>::dft(const auto& Point) 
{
    std::complex<double> FourierTransform = 0;
    std::complex<double> im2pi = im*2*pi;
    for(int i=0; i<Mesh.size(); i++){
        double DotProduct = 0;
        for(int ix=0; ix<3; ix++){
            DotProduct += Point[ix]*Mesh[i][ix];
        }
        FourierTransform += std::exp(im2pi*Dotproduct)*InputArray[i];
    }
    return FourierTransform;
}


