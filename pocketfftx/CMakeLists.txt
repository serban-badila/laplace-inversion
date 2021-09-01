# this is copied into the fetched project source directory.

add_library(pocketfftx INTERFACE)
target_include_directories(pocketfftx INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})