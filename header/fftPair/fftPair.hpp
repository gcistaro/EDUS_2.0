#include <vector>
#include <memory>
#include <complex>
#include "mdContainers/mdContainers.hpp"
#include <fftw3.h>


//NOTE WARNING!!!
//By construction, this class does not own the object InputArray (it does not destruct it) but it owns OutputArray.

//Dimension of InputArray and OutputArray: [howmany x TotalSize of FFt]
class FourierTransform
{
    private:
        mdarray<std::complex<double>, 2> InputArray;
        mdarray<std::complex<double>, 2> OutputArray;
        std::vector< std::vector<double> > Mesh;
        std::vector< std::vector<double> > Mesh_fft;
        int dim;
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
            dim = Dimensions__.size();
            Mesh.resize(dim);
            Mesh_fft.resize(dim);
            InputArray = mdarray<std::complex<double>, 2>(&Data_[0], Data_.get_Size());
            TotalSize = InputArray.get_Size()[1];
            idist = TotalSize;
            odist = TotalSize;
            howmany = Data_.get_Size(0);
            Dimensions = Dimensions__;
            OutputArray = mdarray<std::complex<double>, 2>(InputArray.get_Size());
                
                //Mesh[i].resize(Dimensions[i]);
                //Mesh_fft[i].resize(Dimensions[i]);
                ////by default, we sample always [0,Npoints) in direct space, i.e. [0,1) in fourier space
                //for(int ix=0; ix<Dimensions[i]; ix++){
                //    Mesh[i][ix] = double(ix);
                //    Mesh_fft[i][ix] = double(ix)/double(Dimensions[i]);
                //}
        };

        //FourierTransform(const mdarray<std::complex<double>, 2>& Input__, const mdarray<std::complex<double>, 2>& Output__)
        //{
        //    initialize(Input__, Output__);
        //}

        //void initialize(const mdarray<std::complex<double>, 2>& Input__, const mdarray<std::complex<double>, 2>& Output__)
	    //{
        //    InputArray = std::make_shared<mdarray<std::complex<double>, 2>>(Input__);
        //    OutputArray = std::make_shared<mdarray<std::complex<double>, 2>>(Output__);
        //    howmany = InputArray->get_Size(0);
        //    TotalSize = InputArray->get_Size(1);
        //    dim = Dimensions__.size();
        //    idist = TotalSize;
        //    odist = TotalSize;
        //}

        void initialize(mdarray<std::complex<double>, 2>& Input__, const std::vector<std::vector<double>>& Mesh__)
        {
            Mesh = Mesh__;
            dim = Mesh[0].size();
	        howmany = Input__.get_Size(0);
	        TotalSize = Input__.get_Size(1);
            InputArray = mdarray<std::complex<double>, 2>(&Input__[0], Input__.get_Size());
        }

        void fft(const int& sign);
        mdarray<std::complex<double>, 1> dft(const std::vector<double>& Point, const int& sign); 
        mdarray<std::complex<double>, 2>& dft(const std::vector<std::vector<double>>& ArrayOfPoints, const int& sign);
        
        const mdarray<std::complex<double>, 2>& get_OutputArray() const
        {
            return (OutputArray);
        }
};


void FourierTransform::fft(const int& sign)
{
    //TODO: be sure Data is allocated till the end
    MyPlan = fftw_plan_many_dft(dim, //Dimension of the fft
	                       &Dimensions[0],   
                           howmany,       
                           reinterpret_cast<fftw_complex*>(&InputArray[0]),                    //row-major ordered 
                           inembed,
                           istride,
                           idist,
                           reinterpret_cast<fftw_complex*>(&OutputArray[0]),                   //row-major ordered 
                           onembed,
                           ostride,
                           odist,
                           sign,
			               FFTW_ESTIMATE);
    fftw_execute(MyPlan);
    fftw_destroy_plan(MyPlan);
}


mdarray<std::complex<double>, 1> FourierTransform::dft(const std::vector<double>& Point, const int& sign) 
{
    assert(Point.size() == dim);
    mdarray<std::complex<double>, 1> FT({InputArray.get_Size()[0]});
    FT.fill(std::complex<double>(0.));

    //std::complex<double> FourierTransform = 0;
    static std::complex<double> im2pi = im*2.*pi;
    for(int h=0; h<howmany; h++){
        for(int i=0; i<Mesh.size(); i++){
            double DotProduct = 0;
            for(int ix=0; ix<dim; ix++){
                DotProduct += Point[ix]*Mesh[i][ix];
            }
            FT(h) += std::exp(double(sign)*im2pi*DotProduct)*InputArray(h,i);
	    std::cout<< "DotProduct" << DotProduct << "std::exp(double(sign)*im2pi*DotProduct) " << std::exp(double(sign)*im2pi*DotProduct) << std::endl;//std::cout << "h "<< h <<" howmany " << howmany << "Mesh.size() "<< Mesh.size()<< " i "<<i <<  " FT(h) " << FT(h) << "(*inputarray(h,i) " << (*InputArray)(h,i)<< std::endl;
	}
    }
    return FT;
}

mdarray<std::complex<double>, 2>& FourierTransform::dft(const std::vector<std::vector<double>>& ArrayOfPoints, const int& sign)
{
    OutputArray = mdarray<std::complex<double>, 2>({InputArray.get_Size()[0], ArrayOfPoints.size()});
    for(int ip=0; ip<ArrayOfPoints.size(); ++ip){
        auto FT = dft(ArrayOfPoints[ip], sign);
        //copy to OutputArray
        for(int h=0; h<howmany; ++h){
            OutputArray(h,ip) = FT(h);
            std::cout << "p: " << ip <<  "h " << h << "(*OutputArray)(h,ip) " << (OutputArray)(h,ip) << std::endl;
	}
    }
    return OutputArray;
}


