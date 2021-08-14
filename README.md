# laplace-inversion
Numerical inversion of Laplace transforms using Den Iseger's algorithm [1].

# Usage

See the [demo.py](py_laplace_inversion/demo.py)

# Dev guide

Not yet tested on platforms other than Linux for now.

## Build with cmake 
To build this project (cmake) and run the cxx tests

```
cmake -S . -B cxx-build && cmake --build ./cxx-build && cd cxx-build && ctest
```

from the project root. This will build all targets including the python extension as a shared object.


## Build the python wheels and run the pytests

In a clean python virtual environment: 

```
pip install -rrequirements-dev.txt
python setup.py bdist_wheel
pip install ./dist/*
```

Besides building (cmake) the shared object, `setup.py` will also copy it into `py_laplace_inversion` and this extension can be imported as any other python module. 
Now you can run the python tests for the generated bindings

```
python -m pytest python_tests/
```


[1] https://doi.org/10.1017/S0269964806060013