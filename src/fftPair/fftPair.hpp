#ifndef FOURIERTRANSFORM_HPP
#define FOURIERTRANSFORM_HPP

#include <vector>
#include <memory>
#include <complex>
#include <fftw3.h>

#include "mdContainers/mdContainers.hpp"
#include "Constants.hpp"
//NOTE WARNING!!!
//dft: By construction, this class does not own the objects Array_x and Array_k

//Dimension of Array_x and Array_k: [howmany x TotalSize of FFt]
class FourierTransform
{
    private:
        mdarray<std::complex<double>, 2>* Array_x;
        mdarray<std::complex<double>, 2>* Array_k;
        std::vector< std::vector<double> > Mesh;
        std::vector< std::vector<double> > Mesh_fft;
        int dim;
        std::vector<int> Dimensions;
        int TotalSize;
        double SqrtTotalSize;
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

        FourierTransform(mdarray<std::complex<double>, 2>& Array_x__, 
                         mdarray<std::complex<double>, 2>& Array_k__, 
                         const std::vector<int>& Dimensions__);

        void initialize(mdarray<std::complex<double>, 2>& Array_x__, 
                        mdarray<std::complex<double>, 2>& Array_k__, 
                        const std::vector<int>& Dimensions__);
        void initialize(mdarray<std::complex<double>, 2>& Array_x__, 
                        const std::vector<std::vector<double>>& Mesh__);
        void fft(const int& sign);
        mdarray<std::complex<double>, 1> dft(const std::vector<double>& Point, const int& sign); 
        mdarray<std::complex<double>, 2> dft(const std::vector<std::vector<double>>& ArrayOfPoints, const int& sign);
        inline const mdarray<std::complex<double>, 2>& get_Array_k() const { return (*Array_k);};
};

#endif