project(siprosV4)
cmake_minimum_required(VERSION 3.1)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

message("${CMAKE_BUILD_TYPE}")
message("${CMAKE_CXX_COMPILER}")

set(CMAKE_CXX_STANDARD 17)
find_package(OpenMP REQUIRED)

aux_source_directory(${PROJECT_SOURCE_DIR}/src SRC_LIST)
include_directories(${PROJECT_SOURCE_DIR}/include)

SET(EXECUTABLE_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/bin")

add_subdirectory(openmp)
add_subdirectory(mpi)
add_subdirectory(test)

