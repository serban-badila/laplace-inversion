cmake_minimum_required(VERSION 3.14)

set(CMAKE_VERBOSE_MAKEFILE OFF)

# set the project name
project(LaplaceInversion VERSION 0.1)

# at least c++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)


# include fftw
include(ExternalProject)

MESSAGE(STATUS "Installing fftw3...")
ExternalProject_Add(project_fftw
  URL http://www.fftw.org/fftw-3.3.9.tar.gz
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/fftw
  CONFIGURE_COMMAND ${CMAKE_CURRENT_BINARY_DIR}/fftw/src/project_fftw/configure
  --prefix=${CMAKE_CURRENT_BINARY_DIR}/fftw/install
  INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/fftw/install
  )

add_library(fftw SHARED STATIC IMPORTED)
set(lib_fftw_name ${CMAKE_STATIC_LIBRARY_PREFIX}fftw3${CMAKE_STATIC_LIBRARY_SUFFIX})
set_target_properties(
  fftw 
  PROPERTIES IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/fftw/install/lib/${lib_fftw_name}
)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/fftw/install/include)


add_executable(LaplaceInversion src/laplaceInversion.cpp)
add_dependencies(LaplaceInversion project_fftw)
