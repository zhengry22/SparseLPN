#include"SparseMatrix.h"

using namespace std;
using namespace NTL;
//#define USING_OPENMP

SparseMatrix::SparseMatrix(long row, long col, ZZ q): rows(row), cols(col), q(q) {}
SparseMatrix::SparseMatrix(long row, long col, long k, ZZ q): rows(row), cols(col), sparsity(k), q(q) {}
SparseMatrixCSR::SparseMatrixCSR(long row, long col, ZZ q): SparseMatrix(row, col, q) {
    this->row_ptr.assign(row + 1, 0); 
    this->col_indices.clear();
    this->values.clear();
}
SparseMatrixCSR::SparseMatrixCSR(long row, long col, long k, ZZ q): SparseMatrix(row, col, k, q) {
    this->row_ptr.assign(row + 1, 0); 
    this->col_indices.clear();
    this->values.clear();

    this->col_indices.reserve(row * k);
    this->values.reserve(row * k);
}
long SparseMatrix::getrows() {return this->rows;}
long SparseMatrix::getcols() {return this->cols;}
long SparseMatrix::getsparsity() {return this->sparsity;}
ZZ SparseMatrix::getmod() {return this->q;}

std::unique_ptr<SparseMatrix> SparseMatrixCSRSampler::sample_RDiag(long rows, long cols, long k, const ZZ& q) {
    auto mat = std::make_unique<SparseMatrixCSR>(rows, cols, k, q);

    // mat->row_ptr.assign(rows + 1, 0); 
    // mat->col_indices.clear();
    // mat->values.clear();

    // mat->col_indices.reserve(rows * k);
    // mat->values.reserve(rows * k);

    // 建议使用类成员变量中的 gen 和 seed，这里为了演示保留 local
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<long> dis(0, cols - 1);

    for (long i = 0; i < rows; ++i) {
        // 注意：论文中索引通常从 1 开始，C++ 从 0 开始
        // 这里的 i 代表行索引。论文定义：if i mod (n+1) != 0
        // 我们假设逻辑行号为 i_logic = i + 1
        long i_logic = i + 1; 
        long diagonal_col = i_logic % (cols + 1);
        
        vector<long> current_row_cols;

        // --- 核心修改：处理强制非零元 ---
        if (diagonal_col != 0) {
            // 论文要求第 (diagonal_col - 1) 列必须非零 (调整回 0-indexed)
            current_row_cols.push_back(diagonal_col - 1);
        }

        // --- 采样剩余的 k - (是否已有对角元) 个位置 ---
        while (current_row_cols.size() < (size_t)k) {
            long col = dis(gen);
            // 确保不与强制位置重复，也不与已采样的重复
            if (find(current_row_cols.begin(), current_row_cols.end(), col) == current_row_cols.end()) {
                current_row_cols.push_back(col);
            }
        }
        
        // CSR 格式必须对列索引排序
        sort(current_row_cols.begin(), current_row_cols.end());

        // 填充数据
        for (long col : current_row_cols) {
            mat->col_indices.push_back(col);
            
            ZZ rand_val;
            // 论文要求特定的位置必须是 non-zero
            do {
                rand_val = RandomBnd(q - 1) + 1;
            } while (rand_val == 0); // 确保是非零元
            
            mat->values.push_back(rand_val);
        }

        mat->row_ptr[i + 1] = mat->col_indices.size();
    }

    return mat;
}

vec_ZZ SparseMatrixCSR::operator*(const vec_ZZ& x) const {
    vec_ZZ y;
    y.SetLength(rows);
#ifdef USING_OPENMP
    #pragma omp parallel for
#endif
    for (long i = 0; i < rows; ++i) {
        clear(y[i]);
        for (long j = row_ptr[i]; j < row_ptr[i+1]; j++) {
            // y[i] += values[j] * x[col_indices[j]]
            ZZ tmp;
            mul(tmp, values[j], x[col_indices[j]]);
            add(y[i], y[i], tmp);
        }
        if (this->q != 0) y[i] %= this->q; 
    }
    return y;
} // 这是最基础版本的 CSR 矩阵向量乘法，后面或许可以优化

std::unique_ptr<SparseMatrix> SparseMatrixCSR::operator*(const ZZ& k) const {
    // 1. 在并行区域外创建 C
    auto C = std::make_unique<SparseMatrixCSR>(rows, cols, q);
    
    // 2. 拷贝 CSR 的结构（行指针和列索引），标量乘法不改变这些
    C->row_ptr = this->row_ptr;
    C->col_indices = this->col_indices;
    
    // 3. 必须 resize 分配空间，否则 C->values[i] 会引发段错误
    C->values.resize(this->values.size());

    // 获取数据指针，方便 OpenMP 操作（有些编译器对 vector 迭代器的并行支持更好）
    long nnz = this->values.size();

    // 4. 并行计算数值
#ifdef USING_OPENMP
    #pragma omp parallel for
#endif
    for (long i = 0; i < nnz; ++i) {
        // 计算 C[i] = values[i] * k mod q
        mul(C->values[i], this->values[i], k);
        if (this->q != 0) {
            C->values[i] %= this->q;
        }
    }

    return C;
}

std::unique_ptr<SparseMatrix> SparseMatrixCSR::operator+(const std::unique_ptr<SparseMatrix>& B) const {
    auto b_ptr = dynamic_cast<SparseMatrixCSR*>(B.get());
    
    // 2. 安全检查：如果 B 实际上不是 CSR 格式，转换会返回 nullptr
    if (!b_ptr) {
        throw std::runtime_error("Addition only supported between two CSR matrices");
    }

    // Assert that the 2 matrices has the same dimension
    if (this->rows != b_ptr->getrows() || this->cols != b_ptr->getcols()) 
        throw runtime_error("Dimension mismatch");

    auto C = std::make_unique<SparseMatrixCSR>(rows, cols, q);
    C->row_ptr[0] = 0;

    for (long i = 0; i < this->rows; ++i) {
        long pA = this->row_ptr[i];
        long endA = this->row_ptr[i+1];
        long pB = b_ptr->row_ptr[i];
        long endB = b_ptr->row_ptr[i+1];

        // 归并单行
        while (pA < endA || pB < endB) {
            if (pA < endA && (pB >= endB || this->col_indices[pA] < b_ptr->col_indices[pB])) {
                // 只有 A 有值或 A 的列更靠前
                C->col_indices.push_back(this->col_indices[pA]);
                C->values.push_back(this->values[pA] % q);
                pA++;
            } 
            else if (pB < endB && (pA >= endA || b_ptr->col_indices[pB] < this->col_indices[pA])) {
                // 只有 B 有值
                C->col_indices.push_back(b_ptr->col_indices[pB]);
                C->values.push_back(b_ptr->values[pB] % q);
                pB++;
            } 
            else {
                // 两者列索引相同，相加
                ZZ sum;
                add(sum, this->values[pA], b_ptr->values[pB]);
                if (sum != 0) { // 只有非零才存储
                    C->col_indices.push_back(this->col_indices[pA]);
                    C->values.push_back(sum % this->q);
                }
                pA++; pB++;
            }
        }
        C->row_ptr[i+1] = C->col_indices.size();
    }
    return C; 
}

std::unique_ptr<SparseMatrix> SparseMatrixCSR::operator*(const std::unique_ptr<SparseMatrix>& B) const {
    auto b_ptr = dynamic_cast<SparseMatrixCSR*>(B.get());
    
    // 2. 安全检查：如果 B 实际上不是 CSR 格式，转换会返回 nullptr
    if (!b_ptr) {
        throw std::runtime_error("Addition only supported between two CSR matrices");
    }

    if (this->cols != b_ptr->getrows()) throw runtime_error("Dimension mismatch");

    auto C = std::make_unique<SparseMatrixCSR>(rows, cols, q);
    vector<long> C_col_indices;
    vector<ZZ> C_values;
    C->row_ptr[0] = 0;

    // SPA (Sparse Accumulator) 辅助结构
    // 用于在处理每一行时进行高效累加
    vector<ZZ> spa_values(B->getcols());
    vector<bool> spa_occupied(B->getcols(), false);
    vector<long> spa_indices; // 记录当前行有哪些列变成了非零

    for (long i = 0; i < this->rows; ++i) {
        spa_indices.clear();

        // 遍历 A 的第 i 行
        for (long j = this->row_ptr[i]; j < this->row_ptr[i+1]; ++j) {
            long a_col = this->col_indices[j];
            const ZZ& a_val = this->values[j];

            // 累加 B 的第 a_col 行
            for (long k = b_ptr->row_ptr[a_col]; k < b_ptr->row_ptr[a_col+1]; ++k) {
                long b_col = b_ptr->col_indices[k];
                const ZZ& b_val = b_ptr->values[k];

                if (!spa_occupied[b_col]) {
                    spa_occupied[b_col] = true;
                    spa_indices.push_back(b_col);
                    clear(spa_values[b_col]);
                }

                ZZ prod;
                mul(prod, a_val, b_val);
                add(spa_values[b_col], spa_values[b_col], prod);
            }
        }

        // 将 SPA 中的结果压缩进 C
        sort(spa_indices.begin(), spa_indices.end()); // 保持 CSR 列有序
        for (long col : spa_indices) {
            if (spa_values[col] != 0) {
                C_col_indices.push_back(col);
                C_values.push_back(spa_values[col] % this->q);
                spa_occupied[col] = false; // 重置状态供下一行使用
            }
        }
        C->row_ptr[i+1] = C_col_indices.size();
    }

    C->col_indices = C_col_indices;
    C->values = C_values;
    return C;
}

NTL::mat_ZZ ToDense(SparseMatrixCSR& sparse) {
    NTL::mat_ZZ dense;
    
    // 1. 初始化维度
    dense.SetDims(sparse.rows, sparse.cols);

    // 2. 遍历行指针
    // 注意：即便 values 为空，dense 也会保持全零状态
    for (long i = 0; i < sparse.rows; ++i) {
        // row_ptr[i] 是当前行在 values 数组中的起始位置
        // row_ptr[i+1] 是结束位置（不含）
        for (long j = sparse.row_ptr[i]; j < sparse.row_ptr[i + 1]; ++j) {
            long col = sparse.col_indices[j];
            const NTL::ZZ& val = sparse.values[j];
            
            dense[i][col] = val;
        }
    }

    return dense;
}