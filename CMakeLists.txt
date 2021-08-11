cmake_minimum_required(VERSION 3.14)

set(CMAKE_VERBOSE_MAKEFILE OFF)

project(LaplaceInversion VERSION 0.1)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall --std=c++11 -O3 -fPIC")


# build the FFTW dependency
include(ExternalProject)
ExternalProject_Add(project_fftw
  URL http://www.fftw.org/fftw-3.3.9.tar.gz
  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/fftw
  CONFIGURE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/fftw/src/project_fftw/configure --with-pic
  # The option above will force the configure shell to add the PIC flag (CMake has no control here), 
  # this is needed to be able to link against the generated pybind module later on
  --prefix=${CMAKE_CURRENT_SOURCE_DIR}/fftw/install
  INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/fftw/install
  )

add_library(fftw SHARED STATIC IMPORTED)
set(lib_fftw_name ${CMAKE_STATIC_LIBRARY_PREFIX}fftw3${CMAKE_STATIC_LIBRARY_SUFFIX})
set_target_properties(
  fftw 
  PROPERTIES IMPORTED_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/fftw/install/lib/${lib_fftw_name}
)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/fftw/install/include)

##### Add GTest #####

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(
  laplace_test
  tests/tests.cpp
)

##### Generating the Gtest and pybind11 module targets #####

add_subdirectory(src)

add_dependencies(laplace_test project_fftw)
target_link_libraries(
  laplace_test
  gtest_main
  laplace-inversion
)
target_include_directories(laplace_test PUBLIC "${PROJECT_SOURCE_DIR}/src")
include(GoogleTest)
gtest_discover_tests(laplace_test)


add_subdirectory(pybind11)
add_subdirectory(python_bindings)
