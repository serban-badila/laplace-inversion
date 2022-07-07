![CMake Build](https://github.com/serban-badila/laplace-inversion/actions/workflows/cmake.yml/badge.svg)
![Py Build](https://github.com/serban-badila/laplace-inversion/actions/workflows/pytest.yml/badge.svg)

git@github.com:serban-badila/laplace-inversion.git

# laplace-inversion
Numerical inversion of Laplace transforms using Den Iseger's algorithm [1].

Depends on [pocketfft](https://gitlab.mpcdf.mpg.de/mtr/pocketfft/-/tree/cpp) and [pybind11](https://github.com/pybind/pybind11).


# Usage
Requires a c++11 capable compiler present on the target system and cmake >=3.14 (still have to test different compilers).

See the [demo.py](py_laplace_inversion/demo.py)

# Dev guide

Not tested on platforms other than Linux yet.

## Build with cmake 
To build this project (cmake) and run the cxx tests

```
cmake -S . -B cxx-build/debug -DCMAKE_BUILD_TYPE=Debug && cmake --build ./cxx-build/debug --config Debug && cd cxx-build/debug && ctest
```

from the project root. This will build all targets including the python extension as a shared object.


## Build the python wheels and run the pytests

In a clean python virtual environment: 

```
pip install -rrequirements-dev.txt
python setup.py bdist_wheel
pip install ./dist/*
```

Besides building (cmake) the shared object in a separate build tree (_skbuild), `setup.py` will also copy it into `py_laplace_inversion` and this extension can be imported as any other python module. 
Now you can run the python tests for the generated bindings

```
python -m pytest python_tests/
```

# TODO
 - multi-dimensional inversion (see also [2])
 - ensure you can build this library on any platform and start distributing the binaries 
 - try to implement the general inversion algorithm (for handling functions with arbitrary discontinuities and singularities)

 [1] https://doi.org/10.1017/S0269964806060013

 [2] http://dx.doi.org/10.2139/ssrn.1355449 