cmake_minimum_required(VERSION 3.10)

# 设置项目名称和使用的 C++ 标准
project(Doublepack_Project LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# 设置编译器选项
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O2 -pthread -march=native")

# 添加 NTL, GMP 和数学库的链接
find_library(NTL_LIBRARY ntl REQUIRED)
find_library(GMP_LIBRARY gmp REQUIRED)
find_library(M_LIBRARY m REQUIRED)

# 定义各个源文件
set(SOURCES
    PSS/pss.cpp
    GR/gr.cpp
    RMFE/RMFE.cpp
    RMFE/HenselLift.cpp
    RMFE/NM_RMFE.cpp
    RMFE/utils.cpp
    RMFE/Inverse.cpp
)

# 添加测试文件
add_executable(pss test/pss_test.cpp ${SOURCES})
add_executable(grmult test/GRmult_test.cpp ${SOURCES})
add_executable(rmfe test/RMFE_test.cpp ${SOURCES})

# 链接外部库
target_link_libraries(pss ${NTL_LIBRARY} ${GMP_LIBRARY} ${M_LIBRARY})
target_link_libraries(grmult ${NTL_LIBRARY} ${GMP_LIBRARY} ${M_LIBRARY})
target_link_libraries(rmfe ${NTL_LIBRARY} ${GMP_LIBRARY} ${M_LIBRARY})
