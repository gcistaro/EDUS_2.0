
#include <cuda_runtime.h>
#include <cuComplex.h>
#include <complex>

__global__
void axpby_kernel(
    cuDoubleComplex* __restrict__ Output,
    const cuDoubleComplex* a,
    const cuDoubleComplex* __restrict__ x,
    const cuDoubleComplex* b,
    const cuDoubleComplex* __restrict__ y, 
    int N)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;

    Output[i] = cuCadd(cuCmul(*a, x[i]),
                       cuCmul(*b, y[i]));
}




void axpby_gpu(std::complex<double>* Output, 
               const std::complex<double> a, 
               const std::complex<double>* x, 
               const std::complex<double> b, 
               const std::complex<double>* y,
               int N)
{
    int threads = 256;
    int blocks  = (N+threads-1)/threads;
    axpby_kernel<<<threads, blocks>>>(reinterpret_cast<cuDoubleComplex*>(Output),
                                      reinterpret_cast<const cuDoubleComplex*>(&a),
                                      reinterpret_cast<const cuDoubleComplex*>(x),
                                      reinterpret_cast<const cuDoubleComplex*>(&b),
                                      reinterpret_cast<const cuDoubleComplex*>(y),
                                      N);
}
