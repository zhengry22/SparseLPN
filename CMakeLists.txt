cmake_minimum_required(VERSION 3.10)
project(SHE_Project CXX)

set(CMAKE_CXX_STANDARD 17)

# 1. 查找依赖库
find_library(NTL_LIB ntl REQUIRED)
find_library(GMP_LIB gmp REQUIRED)

# 查找 Eigen3 (这是标准写法)
find_package(Eigen3 REQUIRED)

# 2. 包含目录
include_directories(include)
# 将 Eigen3 的头文件路径加入编译
include_directories(${EIGEN3_INCLUDE_DIR})

# 3. 定义核心库
add_library(she_core
    src/LinearHE.cpp
    src/Matrix.cpp
    src/polynomial.cpp
    src/EncryptionScheme.cpp
)
# 链接 NTL 和 GMP
target_link_libraries(she_core ${NTL_LIB} ${GMP_LIB})

# 4. 编译主程序
add_executable(she_main src/main.cpp)
target_link_libraries(she_main she_core)

# 5. 编译测试程序
add_executable(testLHE tests/testLHE.cpp)
target_link_libraries(testLHE she_core)