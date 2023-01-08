set(CMAKE_C_COMPILER "mpicc")
set(CMAKE_CXX_COMPILER "mpicxx")

add_executable(SiprosEnsembleMPI ${SRC_LIST} ${PROJECT_SOURCE_DIR}/mpi/mpimain.cpp)
target_link_libraries(SiprosEnsembleMPI PUBLIC OpenMP::OpenMP_CXX MStoolkit)
target_compile_options(SiprosEnsembleMPI PUBLIC -ffast-math -march=native -fPIC)
target_compile_definitions(SiprosEnsembleMPI PUBLIC GCC _FILE_OFFSET_BITS=64)