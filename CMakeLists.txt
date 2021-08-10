cmake_minimum_required(VERSION 3.14)

set(CMAKE_VERBOSE_MAKEFILE OFF)

# set the project name
project(LaplaceInversion VERSION 0.1)

# at least c++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall --std=c++11 -O3")


# include fftw
include(ExternalProject)

MESSAGE(STATUS "Installing fftw3...")
ExternalProject_Add(project_fftw
  URL http://www.fftw.org/fftw-3.3.9.tar.gz
  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/fftw
  CONFIGURE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/fftw/src/project_fftw/configure
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


add_executable(LaplaceInversion src/laplaceInversion.cpp)
add_dependencies(LaplaceInversion project_fftw)

target_link_libraries(LaplaceInversion
  fftw
)



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

include_directories(${GTEST_INCLUDE_DIRS} tests)

add_executable(
  laplace_test
  src/laplaceInversion.h
  src/tests.cpp
)

# add_dependencies(laplace_test project_fftw)

target_include_directories(laplace_test PUBLIC "${PROJECT_BINARY_DIR}/src")

target_link_libraries(
  laplace_test
  gtest_main
  fftw
)

include(GoogleTest)
gtest_discover_tests(laplace_test)