add_library(laplace-inversion laplaceInversion.h laplaceInversion.cpp)
target_link_libraries(laplace-inversion pocketfftx)
