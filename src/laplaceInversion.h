#include <vector>
#include <complex>
#include <functional>

#include <fftw3.h>

#define PI 3.141592653589793238462643383279502884
#define OVRSMPL 8

namespace laplaceInversion {

std::vector<double> oneDimensionalInverse(
    std::function<std::complex<double>(std::complex<double>)> f_hat,
    double delta, 
    unsigned int mexp,
    int n
);

} // namespace laplaceInversion