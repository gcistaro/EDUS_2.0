#include <vector>
#include <memory>
#include <complex>
#include "../mdContainers/mdContainers.hpp"
#include <fftw3.h>


//Dimension of InputArray and OutputArray: [howmany x TotalSize of FFt]
class FourierTransform
{
    private:
        std::shared_ptr<mdarray<std::complex<double>, 2>> InputArray;
        std::shared_ptr<mdarray<std::complex<double>, 2>> OutputArray;
        std::vector<std::array<double>, 3> Mesh;
        std::vector<std::array<double>, 3> Mesh_fft;
        std::vector<int> Dimensions;
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
        FourierTransform(){};
        FourierTransform(mdarray<std::complex<double>,2>& Data_, std::vector<int>& Dimensions__)
        {
            Mesh.resize(dim);
            Mesh_fft.resize(dim);
            InputArray = std::make_shared<mdarray<std::complex<double>, 2>>(Data_);
            TotalSize = InputArray->get_Size()[1];
            idist = TotalSize;
            odist = TotalSize;

            Dimensions = Dimensions__;
            OutputArray = std::make_shared<mdarray<std::complex<double>, 2>>(
                                mdarray<std::complex<double>, dim>(InputArray->get_Size()));
                
                //Mesh[i].resize(Dimensions[i]);
                //Mesh_fft[i].resize(Dimensions[i]);
                ////by default, we sample always [0,Npoints) in direct space, i.e. [0,1) in fourier space
                //for(int ix=0; ix<Dimensions[i]; ix++){
                //    Mesh[i][ix] = double(ix);
                //    Mesh_fft[i][ix] = double(ix)/double(Dimensions[i]);
                //}
        };

        FourierTransform(const mdarray<std::complex<double>, 2>& Input__, const mdarray<std::complex<double>, 2>& Output__)
        {
            initialize(Input__, Output__);
        }

        void initialize(const mdarray<std::complex<double>, 2>& Input__, const mdarray<std::complex<double>, 2>& Output__)
	    {
            InputArray = std::make_shared<mdarray<std::complex<double>, 2>>(Input__);
            OutputArray = std::make_shared<mdarray<std::complex<double>, 2>>(Output__);
        }

        FourierTransform(const BlockMatrix<T,space>& BlockMatrix__)
        {
            InputArray = std::make_shared<mdarray<std::complex<double>, 2>>(BlockMatrix__.get_Values());
            Mesh.resize(BlockMatrix__.get_mesh().size());
        }



        void fft(const int& sign);
        mdarray<std::complex<double>, 1> dft(const std::array<double,dim>& Point, const int& sign); 
        mdarray<std::complex<double>, 2>& dft(const std::vector<std::array<double,dim>>& ArrayOfPoints, const int& sign);
        
        const mdarray<std::complex<double>, dim>& get_OutputArray() const
        {
            return (*OutputArray);
        }
};


template<size_t dim>
void FourierTransform<dim>::fft(const int& sign)
{
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
mdarray<std::complex<double>, 1> FourierTransform<dim>::dft(const std::array<double,dim>& Point, const int& sign) 
{
    mdarray<std::complex<double>, 1> FT(InputArray->get_Size());
    FT.fill(std::complex<double>(0.));

    //std::complex<double> FourierTransform = 0;
    static std::complex<double> im2pi = im*2.*pi;
    for(int h=0; h<howmany; h++){
        for(int i=0; i<Mesh.size(); i++){
            double DotProduct = 0;
            for(int ix=0; ix<dim; ix++){
                DotProduct += Point[ix]*Mesh[i][ix];
            }
            FT(h) += std::exp(double(sign)*im2pi*DotProduct)*(*InputArray)(h,i);
        }
    }
    return FT;
}

template<size_t dim>
mdarray<std::complex<double>, 2>& FourierTransform<dim>::dft(const std::vector<std::array<double,dim>>& ArrayOfPoints, const int& sign)
{
    OutputArray = std::make_shared<mdarray<std::complex<double>, 2>>(mdarray<std::complex<double>, dim>({InputArray->get_Size()[0], ArrayOfPoints.size()}));
    for(int ip=0; ip<ArrayOfPoints.size(); ++ip){
        auto FT = dft(ArrayOfPoints[ip], sign);
        //copy to OutputArray
        for(int h=0; h<howmany; ++h){
            (*OutputArray)(h,ip) = FT(ip);
        }
    }
    return *OutputArray;
}


