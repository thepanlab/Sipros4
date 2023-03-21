set(CMAKE_C_COMPILER "mpicc")
set(CMAKE_CXX_COMPILER "mpicxx")

add_executable(SiprosV4mpi ${SRC_LIST} ${PROJECT_SOURCE_DIR}/mpi/mpimain.cpp)
target_link_libraries(SiprosV4mpi PUBLIC OpenMP::OpenMP_CXX)
# target_compile_options(SiprosV4mpi PUBLIC -ffast-math -march=native)
target_compile_options(SiprosV4mpi PUBLIC -ffast-math)