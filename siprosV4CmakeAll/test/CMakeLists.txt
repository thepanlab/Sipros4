
add_executable(SiprosV4test ${SRC_LIST} ${PROJECT_SOURCE_DIR}/test/main.cpp)

target_link_libraries(SiprosV4test PUBLIC OpenMP::OpenMP_CXX)
target_compile_options(SiprosV4test PUBLIC -ffast-math -march=native)