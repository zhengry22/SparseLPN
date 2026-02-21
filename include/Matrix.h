// This file is an implementation of the sparse matrix
#include<iostream>
#include<vector>

class MySparseMatrix {
private:
    std::vector<long> values;     // 非零值
    std::vector<int> col_indices; // 列索引
    std::vector<int> row_ptr;     // 行偏移
    long q;                       // Modulus
public:

};