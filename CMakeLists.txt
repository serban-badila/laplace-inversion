cmake_minimum_required(VERSION 3.14)

set(CMAKE_VERBOSE_MAKEFILE OFF)

project(LaplaceInversion VERSION 0.1)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall --std=c++11 -O3 -fPIC")


include(FetchContent)

##### Fetch pocketfft and configure it with a minimal cmake file ####

MESSAGE(STATUS "Fetching pocketfft...")
FetchContent_Declare(
  pocketfftx 
  GIT_REPOSITORY https://gitlab.mpcdf.mpg.de/mtr/pocketfft.git
  GIT_TAG ea778e37710c07723435b1be58235996d1d43a5a  # on the "cpp" branch
)

FetchContent_Populate(pocketfftx)
FetchContent_GetProperties(pocketfftx)
file(COPY pocketfftx/CMakeLists.txt DESTINATION ${pocketfftx_SOURCE_DIR})
add_subdirectory(${pocketfftx_SOURCE_DIR})


##### Fetch pybind11 #####

MESSAGE(STATUS "Fetching pybind11...")
FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11
    GIT_TAG v2.7.1
)
FetchContent_MakeAvailable(pybind11)

add_subdirectory(src)

if (CMAKE_BUILD_TYPE MATCHES Debug)
  
  ##### Fetch GTest #####
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

  target_link_libraries(
    laplace_test
    gtest_main
    laplace-inversion
    pocketfftx
  )
  target_include_directories(laplace_test PUBLIC "${PROJECT_SOURCE_DIR}/src")
  include(GoogleTest)
  gtest_discover_tests(laplace_test)

endif(CMAKE_BUILD_TYPE MATCHES Debug)

# this requires pybind11 configured 
if(CMAKE_BUILD_TYPE MATCHES Release)
  add_subdirectory(py_laplace_inversion)
endif(CMAKE_BUILD_TYPE MATCHES Release)
