#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>

#include "laplaceInversion.h"

#define DEFAULT_QUAD_SIZE 48


PYBIND11_MODULE(_py_laplace, m) {
    m.doc() = "Compute the inverse of a Laplace(-Stieltjes) Transform using den Iseger's algorithm.";
    
    m.def(
        "laplace_inversion", 
        &laplaceInversion::oneDimensionalInverse,
        pybind11::arg("f"),
        pybind11::arg("delta"),
        pybind11::arg("mexp"),
        pybind11::arg("n") = DEFAULT_QUAD_SIZE, 
        R"pbdoc(
            Compute the inverse of a one-dimensional Laplace(-Stieltjes) Transform using den Iseger's algorithm [1].
        
        Args:
            f (Callable[[complex], complex]): the Laplace tansform to be inverted.
            mexp (int): the exponent defining the number of evaluation points (2^mexp) for the inverse where precision is guaranteed; 
                The transform evaluation uses an oversampling factor, 8 * 2^mexp.
            delta (float): mesh width.
            n (int): The size of the Gauss quadrature; Choose between 16, 32 and 48; Defaults to 48 (recommended for higher precision). 
        
        Returns:
                A list containing the inverse function evaluated at {0, delta, 2*delta, ...}.

        [1] https://papers.ssrn.com/sol3/papers.cfm?abstract_id=1013507
        )pbdoc"
    );
}
