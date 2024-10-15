cmake_minimum_required(VERSION 3.14)

# 设置项目名称和使用的 C++ 标准
project(Doublepack_Project DESCRIPTION "Double-Packing")

# 设置编译器选项
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -pthread -std=gnu++17")

set(FI "fi.x")
set(OURS "ours.x")

# 添加 NTL, GMP 和数学库的链接
find_library(NTL_LIBRARY ntl REQUIRED)
find_library(GMP_LIBRARY gmp REQUIRED)
find_library(M_LIBRARY m REQUIRED)

# 定义各个源文件
set(DP_SOURCE_FILES
    PSS/pss.cpp
    GR/gr.cpp
    RMFE/RMFE.cpp
    RMFE/HenselLift.cpp
    RMFE/NM_RMFE.cpp
    RMFE/utils.cpp
    RMFE/Inverse.cpp

    src/dp/gate.cc
    src/dp/mult_gate.cc
    src/dp/input_gate.cc
    src/dp/output_gate.cc
  
    src/dp/circuits/building.cc
    src/dp/circuits/cleartext.cc
    src/dp/circuits/generic.cc
    src/dp/circuits/online.cc
    src/dp/circuits/offline.cc  
  
    src/dp/correlator.cc  
    src/dp/fi_prep.cc  
)

set(TEST_SOURCE_FILES
    test/main.cc
    test/pss_test.cc
    test/GRmult_test.cc
    test/RMFE_test.cc
    test/fi_prep_test.cc
    test/input_test.cc
    test/output_test.cc
    test/mult_test.cc
    test/channel_specialized_test.cc  # Add this line
)

add_subdirectory( secure-computation-library )

set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -O0" )

set(TEST_EXECUTABLE "tests.x")
find_package( Catch2 REQUIRED )
include( CTest )
include( Catch )

## Executable for tests
add_executable( ${TEST_EXECUTABLE} ${DP_SOURCE_FILES} ${TEST_SOURCE_FILES} )
target_include_directories(${TEST_EXECUTABLE} PUBLIC
    "${CMAKE_SOURCE_DIR}/secure-computation-library/include"
    "${CMAKE_SOURCE_DIR}/RMFE"
    "${CMAKE_SOURCE_DIR}/GR"
    "${CMAKE_SOURCE_DIR}/PSS"
    "${CMAKE_SOURCE_DIR}/src")
set( libscl "${CMAKE_SOURCE_DIR}/secure-computation-library/build/libscl.so.0.3" )
target_link_libraries( ${TEST_EXECUTABLE} Catch2::Catch2 pthread scl ${NTL_LIBRARY} ${GMP_LIBRARY} ${M_LIBRARY} )

#set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -O2")
  ## Executable of our protocol
#  add_executable(${FI} experiment/fi.cc ${DP_SOURCE_FILES})
#  target_include_directories(${FI} PUBLIC
#    "${CMAKE_SOURCE_DIR}/secure-computation-library/include"
#    "${CMAKE_SOURCE_DIR}/RMFE"
#    "${CMAKE_SOURCE_DIR}/GR"
#    "${CMAKE_SOURCE_DIR}/PSS"
#    "${CMAKE_SOURCE_DIR}/src")
#  target_link_libraries(${FI} pthread scl ${NTL_LIBRARY} ${GMP_LIBRARY} ${M_LIBRARY})

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
