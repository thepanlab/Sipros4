# set(CMAKE_C_COMPILER "gcc-11")
# set(CMAKE_CXX_COMPILER "g++-11")

add_executable(SiprosEnsembleOMP ${SRC_LIST} ${PROJECT_SOURCE_DIR}/openmp/main.cpp)

# target_compile_options(SiprosEnsembleOMP PUBLIC -ffast-math -march=native)
target_compile_options(SiprosEnsembleOMP PUBLIC -ffast-math)

target_link_libraries(SiprosEnsembleOMP PUBLIC OpenMP::OpenMP_CXX MStoolkit gomp dl -static)

target_compile_definitions(SiprosEnsembleOMP PUBLIC GCC _FILE_OFFSET_BITS=64)