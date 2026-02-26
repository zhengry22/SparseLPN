// This file is for testing whether we correctly implemented the Sparse Matrix.
#pragma once
#include <NTL/matrix.h>
#include <NTL/vec_ZZ.h>
#include <memory>
#include <NTL/mat_ZZ.h>
#include "SparseMatrix.h" // 引入你的 SparseMatrixCSR 定义

using namespace NTL;

class DenseMatrixZZ {
private:
    mat_ZZ data;
    long rows;
    long cols;
    ZZ q; // 模数

public:
    DenseMatrixZZ(long r, long c, const ZZ& modulus = ZZ::zero());
    
    // 从已有的 NTL 矩阵构造
    DenseMatrixZZ(const mat_ZZ& m, const ZZ& modulus = ZZ::zero());

    // 基础算术
    DenseMatrixZZ operator+(const DenseMatrixZZ& other) const;
    DenseMatrixZZ operator*(const DenseMatrixZZ& other) const;
    vec_ZZ operator*(const vec_ZZ& v) const; // 矩阵-向量乘法
    DenseMatrixZZ operator*(const ZZ& scalar) const; // 数乘

    // 核心接口：转化为 CSR 稀疏矩阵
    std::unique_ptr<SparseMatrixCSR> toCSR() const;

    // 辅助功能
    void set(long r, long c, const ZZ& val);
    const ZZ& get(long r, long c) const;
    long getRows() const { return rows; }
    long getCols() const { return cols; }
    void print();
    
    // 判定相等 (用于测试)
    bool operator==(const DenseMatrixZZ& other) const {
        return (data == other.data && q == other.q);
    }
};