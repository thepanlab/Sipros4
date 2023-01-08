
add_executable(SiprosV3test ${SRC_LIST} ${PROJECT_SOURCE_DIR}/test/main.cpp)

target_link_libraries(SiprosV3test PUBLIC OpenMP::OpenMP_CXX)
target_compile_options(SiprosV3test PUBLIC -ffast-math -march=native)