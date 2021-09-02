pybind11_add_module(_py_laplace laplace_bindings.cpp)

target_link_libraries(_py_laplace PRIVATE laplace-inversion)
target_include_directories(_py_laplace PUBLIC "${PROJECT_SOURCE_DIR}/src")
install(TARGETS _py_laplace CONFIGURATIONS Release LIBRARY DESTINATION py_laplace_inversion)
