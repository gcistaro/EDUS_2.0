#include <vector>
#include <memory>
#include <fftw.h>

template<typename T>
class fftPair
{
    private:
        T* Data    = nullptr;
        T* fftData = nullptr;
    	std::vector<int> Dimensions;
        int Dimension;
        fftw_plan fftw_plan_dft;
    public:
        fftPair(){};
        fftPair(const auto& Dimensions_, const T* Data_)
        {
            Dimensions = Dimensions_;
            Data = Data_;
            Dimension = Dimensions.size();
            
        };
        void set_Data(const T& Data_) {Data = Data_;};
        void set_Dimension(const auto& Dimensions_) {Dimensions = Dimensions_;};
        
	    void fft(const int& sign);

}


template<typename T>
void fftPair<T>::fft(const int& sign)
{
    //TODO: be sure Data is allocated till the end
    fftw_plan_dft = fftw_create_plan(int(Dimensions.size()), //Dimension of the fft
	                                &Dimensions[0],          //array with size of the array over each dimension
                                    Data,                    //row-major ordered 
			                        fftData,
                                    sign,
			                        FFTW_ESTIMATE);
    fftw_execute(fftw_plan_dft);
}


