#include <Eigen/Sparse>
#include <iostream>

typedef __int128 int128; // GCC/Clang 支持的 128 位整数

int main() {
    // 定义稀疏矩阵：值类型为 int128，索引类型为默认的 int
    Eigen::SparseMatrix<int128> mat(100, 100);

    // 填充数据
    mat.insert(0, 0) = (int128)1 << 100; // 一个巨大的 128 位数

    // 运算
    Eigen::SparseMatrix<int128> res = mat * mat;

    return 0;
}