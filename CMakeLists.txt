cmake_minimum_required(VERSION 3.10)
project(SHE_Project CXX)

# ----------------------------------------------------------------              
# 1. 编译模式与全局设置
# ----------------------------------------------------------------
# 设置默认编译模式为 Debug，方便调试。发布性能时改用 -DCMAKE_BUILD_TYPE=Release
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# ----------------------------------------------------------------
# 2. 编译选项 (根据模式自动切换优化等级)
# ----------------------------------------------------------------
if(MSVC)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2 /arch:AVX2")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /Zi /Od")
else()
    # Release: 最高优化级别 + 指令集优化
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native -DNDEBUG")
    
    # Debug: 关闭所有优化，保留完整符号表，开启所有警告
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0 -Wall -Wextra")
endif()

# ----------------------------------------------------------------
# 3. 依赖库查找
# ----------------------------------------------------------------

# 查找 GMP
find_path(GMP_INCLUDE_DIR gmp.h)
find_library(GMP_LIB gmp)
if(NOT GMP_LIB)
    message(FATAL_ERROR "GMP library not found!")
endif()

# 查找 NTL
find_path(NTL_INCLUDE_DIR NTL/ZZ.h)
find_library(NTL_LIB ntl)
if(NOT NTL_LIB)
    message(FATAL_ERROR "NTL library not found!")
endif()

# 查找 Eigen3
find_package(Eigen3 REQUIRED)

# [可选] 如果你以后决定用 OpenMP，取消下面两行的注释
# find_package(OpenMP REQUIRED)
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

# ----------------------------------------------------------------
# 4. 构建核心库 (she_core)
# ----------------------------------------------------------------
file(GLOB CORE_SOURCES "src/*.cpp")

add_library(she_core STATIC ${CORE_SOURCES})

# 使用 target 级别的配置，自动处理依赖传递
target_include_directories(she_core PUBLIC 
    include 
    ${NTL_INCLUDE_DIR} 
    ${GMP_INCLUDE_DIR}
    ${EIGEN3_INCLUDE_DIR}
)

target_link_libraries(she_core PUBLIC 
    ${NTL_LIB} 
    ${GMP_LIB} 
    Eigen3::Eigen
)

# ----------------------------------------------------------------
# 5. 构建可执行文件
# ----------------------------------------------------------------

# 主程序
add_executable(she_main src/main.cpp)
target_link_libraries(she_main she_core)

# 矩阵功能测试程序
add_executable(testSparseMatrix tests/testSparseMatrix.cpp)
target_link_libraries(testSparseMatrix she_core)

# LHE 测试程序
add_executable(testLHE tests/testLHE.cpp)
target_link_libraries(testLHE she_core)

# ----------------------------------------------------------------
# 6. 安装/输出设置 (可选)
# ----------------------------------------------------------------
# 将二进制文件统一输出到 bin 目录
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)