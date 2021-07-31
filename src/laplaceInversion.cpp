#include "laplaceInversion.h"

#include <memory>
#include <iostream>
#include <complex> // include this before the FFT library
#include <cmath>
#include  <vector>

#include <fftw3.h>


/* Test client */

std::complex<double> exp_transform(std::complex<double> s) {
    return 1.0 / (s + .5);
}

double exponential(double x) {
    // This is the theoretical inverse of the above Laplace transform.
    return exp(-.5*x);    
}


int main() {
    int n {48}; // can take only a pre-defined set of values; see the quadrature table
    
    double delta {.1};
    unsigned int mexp {10};
    int m2 = 8 * pow(2, mexp);
    
    std::vector<double> inverse = laplaceInversion::oneDimensionalInverse(&exp_transform, delta, mexp, n);
    
    int index{5};
    std::cout << "Inverse: " << inverse[index] << std::endl;
    std::cout << "Actual value: " << exponential(index * delta) << std::endl;

    return 0;

}