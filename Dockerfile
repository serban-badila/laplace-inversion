FROM python:3.8 as base
    RUN apt update && apt-get install -y cmake

    WORKDIR /app
    ADD requirements-dev.txt .

    RUN pip install -U pip
    RUN pip install -r requirements-dev.txt

FROM base as cxx-build
    WORKDIR /app

    ADD CMakeLists.txt .
    ADD src ./src
    ADD tests ./tests
    ADD pocketfftx ./pocketfftx
    ADD py_laplace_inversion ./py_laplace_inversion
    ADD python_tests ./python_tests
    ADD setup.py .
    ADD VERSION .
    ADD LICENSE .
    ADD README.md .

    RUN cmake -S . -B cxx-build/debug -D CMAKE_BUILD_TYPE=Debug && cmake --build ./cxx-build/debug --config Debug 
    RUN cd cxx-build/debug && ctest

FROM cxx-build as python-build
    RUN python setup.py bdist_wheel
    RUN pip install ./dist/*

FROM python-build as python-tests
    RUN cd .. && python -m pytest app/python_tests
