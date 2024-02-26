#include <vector>
#include <memory>
#include <complex>
#include "mdContainers/mdContainers.hpp"
#include <fftw3.h>


//NOTE WARNING!!!
//dft: By construction, this class does not own the object Array_x (it does not destruct it) but it owns Array_k.
//fft: none of them is owned by the class.

//Dimension of InputArray and Array_k: [howmany x TotalSize of FFt]
class FourierTransform
{
    private:
        mdarray<std::complex<double>, 2> Array_x;
        mdarray<std::complex<double>, 2> Array_k;
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

        FourierTransform(mdarray<std::complex<double>, 2>& Array_x__, mdarray<std::complex<double>, 2>& Array_k__, const std::vector<int>& Dimensions__)
        {
            this->initialize(Array_x__, Array_k__, Dimensions__);
        };

        void initialize(mdarray<std::complex<double>, 2>& Array_x__, mdarray<std::complex<double>, 2>& Array_k__, const std::vector<int>& Dimensions__)
	    {
            Array_x.initialize(&Array_x__(0,0), Array_x__.get_Size());
            Array_k.initialize(&Array_k__(0,0), Array_x__.get_Size());
            howmany = Array_x.get_Size(0);
            TotalSize = Array_x.get_Size(1);
            dim = Dimensions__.size();
            Dimensions = Dimensions__;
            idist = TotalSize;
            odist = TotalSize;
        }

        void initialize(mdarray<std::complex<double>, 2>& Array_x__, const std::vector<std::vector<double>>& Mesh__)
        {
            Mesh = Mesh__;
	        howmany = Array_x__.get_Size(0);
            dim = Mesh[0].size();
	        TotalSize = Array_x__.get_Size(1);
            Array_x = mdarray<std::complex<double>, 2>(&Array_x__[0], Array_x__.get_Size());
        }

        void fft(const int& sign);
        mdarray<std::complex<double>, 1> dft(const std::vector<double>& Point, const int& sign); 
        mdarray<std::complex<double>, 2>& dft(const std::vector<std::vector<double>>& ArrayOfPoints, const int& sign);
        
        const mdarray<std::complex<double>, 2>& get_Array_k() const
        {
            return (Array_k);
        }
};


void FourierTransform::fft(const int& sign)
{
    //TODO: be sure Data is allocated till the end
    auto& input = (sign == +1 ? Array_k : Array_x );//sign = +1 correspond to ifft -> f(x) = sum_n c_n e^{+inx}
    auto& output = (sign == +1 ? Array_x : Array_k); 
    MyPlan = fftw_plan_many_dft(dim, //Dimension of the fft
	                       &Dimensions[0],   
                           howmany,       
                           reinterpret_cast<fftw_complex*>(&input[0]),                    //row-major ordered 
                           inembed,
                           istride,
                           idist,
                           reinterpret_cast<fftw_complex*>(&output[0]),                   //row-major ordered 
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
    mdarray<std::complex<double>, 1> FT({Array_x.get_Size()[0]});
    FT.fill(std::complex<double>(0.));

    //std::complex<double> FourierTransform = 0;
    static std::complex<double> im2pi = im*2.*pi;
    for(int h=0; h<howmany; h++){
        for(int i=0; i<Mesh.size(); i++){
            double DotProduct = 0;
            for(int ix=0; ix<dim; ix++){
                DotProduct += Point[ix]*Mesh[i][ix];
            }
            FT(h) += std::exp(double(sign)*im2pi*DotProduct)*Array_x(h,i);
	    std::cout<< "DotProduct" << DotProduct << "std::exp(double(sign)*im2pi*DotProduct) " << std::exp(double(sign)*im2pi*DotProduct) << std::endl;//std::cout << "h "<< h <<" howmany " << howmany << "Mesh.size() "<< Mesh.size()<< " i "<<i <<  " FT(h) " << FT(h) << "(*Array_x(h,i) " << (*Array_x)(h,i)<< std::endl;
	    }
    }
    return FT;
}

mdarray<std::complex<double>, 2>& FourierTransform::dft(const std::vector<std::vector<double>>& ArrayOfPoints, const int& sign)
{
    Array_k = mdarray<std::complex<double>, 2>({Array_x.get_Size()[0], ArrayOfPoints.size()});
    for(int ip=0; ip<ArrayOfPoints.size(); ++ip){
        auto FT = dft(ArrayOfPoints[ip], sign);
        //copy to Array_k
        for(int h=0; h<howmany; ++h){
            Array_k(h,ip) = FT(h);
            std::cout << "p: " << ip <<  "h " << h << "(*Array_k)(h,ip) " << (Array_k)(h,ip) << std::endl;
	}
    }
    return Array_k;
}


