# set(CMAKE_C_COMPILER "gcc-11")
# set(CMAKE_CXX_COMPILER "g++-11")

add_executable(SiprosEnsembleOMP ${SRC_LIST} ${PROJECT_SOURCE_DIR}/openmp/main.cpp)

# target_compile_options(SiprosEnsembleOMP PUBLIC -ffast-math -march=native)
target_compile_options(SiprosEnsembleOMP PUBLIC -ffast-math)

# target_link_libraries(SiprosEnsembleOMP PUBLIC OpenMP::OpenMP_CXX MStoolkit gomp dl -static)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
message(CMAKE_CXX_FLAGS: "${CMAKE_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
message(CMAKE_EXE_LINKER_FLAGS: "${CMAKE_EXE_LINKER_FLAGS}")
target_link_libraries(SiprosEnsembleOMP PUBLIC MStoolkit)

target_compile_definitions(SiprosEnsembleOMP PUBLIC GCC _FILE_OFFSET_BITS=64)