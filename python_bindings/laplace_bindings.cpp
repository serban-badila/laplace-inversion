#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>

#include "laplaceInversion.h"


namespace py = pybind11;

int foo (int x){
    return 2 * x;
}

PYBIND11_MODULE(py_laplace, m) {

    m.def(
        "laplace_inversion", 
        [](std::function< std::complex<double>(std::complex<double>) > &f, double delta, size_t mexp, int n) {
            return laplaceInversion::oneDimensionalInverse(f, delta, mexp, n);
    }, 
        "Laplace inversion using den Iseger's algorithm"
    );
}
