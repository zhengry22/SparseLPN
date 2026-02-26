#include <NTL/vec_ZZ.h>
#include <NTL/ZZ.h>
#include <vector>
#include <map>
#include <random>
#include <algorithm>
#include <memory>
#include <NTL/mat_ZZ.h>
using namespace std;
using namespace NTL;

class SparseMatrix {
protected:
    long rows, cols; // which corresponds to m, n in the paper
    long sparsity; // sparsity
    ZZ q; // mod
public:
    SparseMatrix(long row, long col, ZZ q);
    SparseMatrix(long row, long col, long k, ZZ q);
    //拷贝构造函数

    // 获取基本信息
    long getrows();
    long getcols();
    long getsparsity();
    ZZ getmod();

    virtual ~SparseMatrix() {}

    // 1. 矩阵-向量乘法: y = A * x
    virtual vec_ZZ operator*(const vec_ZZ& x) const = 0;

    // 2. 矩阵数乘：A' = k * A
    virtual std::unique_ptr<SparseMatrix> operator*(const ZZ& k) const = 0;

    // 3. 矩阵加法：C = A + B
    virtual std::unique_ptr<SparseMatrix> operator+(const std::unique_ptr<SparseMatrix>& B) const = 0;

    // 4. 矩阵乘法：C = A * B
    virtual std::unique_ptr<SparseMatrix> operator*(const std::unique_ptr<SparseMatrix>& B) const = 0;

    // 在稀疏矩阵末尾添加新的一列
    //virtual std::unique_ptr<SparseMatrix> addnewcolumn(const vec_ZZ& b) const = 0;

    // 将稀疏矩阵转化为正常的矩阵形式 (从而进行调试和性能对比)
};

// 这是一个专门用于采样 Diag 矩阵的 Sampler 类
class SparseMatrixSampler {
public: 
    virtual std::unique_ptr<SparseMatrix> sample_RDiag(long n, long m, long k, const ZZ& q) = 0;
};

class SparseMatrixCSR: public SparseMatrix {
public:
    vector<long> row_ptr;       // 长度为 rows + 1
    vector<long> col_indices;   // 长度为 nnz (非零元个数)
    vector<ZZ> values;          // 长度为 nnz

    SparseMatrixCSR(long row, long col, ZZ q);
    SparseMatrixCSR(long row, long col, long k, ZZ q);

    // 1. 矩阵-向量乘法: y = A * x
    vec_ZZ operator*(const vec_ZZ& x) const override;

    // 2. 矩阵数乘: A = k * A
    std::unique_ptr<SparseMatrix> operator*(const ZZ& k) const override;

    // 3. 矩阵加法: C = this + B
    // 假设 A 和 B 维度相同
    std::unique_ptr<SparseMatrix> operator+(const std::unique_ptr<SparseMatrix>& B) const override;

    // 4. 矩阵乘法：return this * B
    std::unique_ptr<SparseMatrix> operator*(const std::unique_ptr<SparseMatrix>& B) const override;

    // 还要添加一个接口，能够将一个向量当中的元素添加到稀疏矩阵的最后一列
    //std::unique_ptr<SparseMatrix> addnewcolumn(const vec_ZZ& b) const override;

    friend NTL::mat_ZZ ToDense(SparseMatrixCSR& sparse);
};

NTL::mat_ZZ ToDense(SparseMatrixCSR& sparse);

// 然后我们给 SparseMatrixCSR 单独写一个生成随机对角阵的类
class SparseMatrixCSRSampler: public SparseMatrixSampler {
public:
    std::unique_ptr<SparseMatrix> sample_RDiag(long n, long m, long k, const ZZ& q) override;
};