pybind11_add_module(py_laplace laplace_bindings.cpp)

target_link_libraries(py_laplace PRIVATE laplace-inversion)
target_include_directories(py_laplace PUBLIC "${PROJECT_SOURCE_DIR}/src")
