add_library(laplace-inversion laplaceInversion.h laplaceInversion.cpp)

add_dependencies(laplace-inversion project_fftw)

target_link_libraries(laplace-inversion fftw)
