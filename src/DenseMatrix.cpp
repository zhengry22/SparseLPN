#include "DenseMatrix.h"

DenseMatrixZZ::DenseMatrixZZ(long r, long c, const ZZ& modulus) 
    : rows(r), cols(c), q(modulus) {
    data.SetDims(r, c);
}

DenseMatrixZZ::DenseMatrixZZ(const mat_ZZ& m, const ZZ& modulus) 
    : data(m), rows(m.NumRows()), cols(m.NumCols()), q(modulus) {}

void DenseMatrixZZ::set(long r, long c, const ZZ& val) {
    if (q > 0) data[r][c] = val % q;
    else data[r][c] = val;
}

const ZZ& DenseMatrixZZ::get(long r, long c) const {
    return data[r][c];
}

// 矩阵加法
DenseMatrixZZ DenseMatrixZZ::operator+(const DenseMatrixZZ& other) const {
    DenseMatrixZZ res(rows, cols, q);
    add(res.data, this->data, other.data);
    if (q > 0) {
        for (long i = 0; i < rows; i++)
            for (long j = 0; j < cols; j++)
                res.data[i][j] %= q;
    }
    return res;
}

// 矩阵乘法 (使用 NTL 自带的高效乘法)
DenseMatrixZZ DenseMatrixZZ::operator*(const DenseMatrixZZ& other) const {
    DenseMatrixZZ res(this->rows, other.cols, q);
    mul(res.data, this->data, other.data);
    if (q > 0) {
        for (long i = 0; i < res.rows; i++)
            for (long j = 0; j < res.cols; j++)
                res.data[i][j] %= q;
    }
    return res;
}

// 数乘
DenseMatrixZZ DenseMatrixZZ::operator*(const ZZ& scalar) const {
    DenseMatrixZZ res(rows, cols, q);
    mul(res.data, this->data, scalar);
    if (q > 0) {
        for (long i = 0; i < rows; i++)
            for (long j = 0; j < cols; j++)
                res.data[i][j] %= q;
    }
    return res;
}

// 核心接口：Dense -> CSR 转换逻辑
std::unique_ptr<SparseMatrixCSR> DenseMatrixZZ::toCSR() const {
    auto csr = std::make_unique<SparseMatrixCSR>(rows, cols, q);
    
    csr->row_ptr.push_back(0);
    for (long i = 0; i < rows; i++) {
        for (long j = 0; j < cols; j++) {
            if (data[i][j] != 0) {
                csr->values.push_back(data[i][j]);
                csr->col_indices.push_back(j);
            }
        }
        csr->row_ptr.push_back(csr->values.size());
    }
    return csr;
}

void DenseMatrixZZ::print() {
    for (long i = 0; i < rows; i++) {
        for (long j = 0; j < rows; j++) {
            cout << this->data[i][j] << " ";
        }
        cout << endl;
    }
}