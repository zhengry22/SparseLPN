#include"SparseMatrix.h"

using namespace std;
using namespace NTL;
#define USING_OPENMP
//#define WITH_ERROR

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
        if (y[i] < 0) y[i] += q;
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
            if (C->values[i] < 0) C->values[i] += q;
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

                ZZ tmp = this->values[pA] % q;
                if (tmp < 0) tmp += q;

                C->values.push_back(tmp);
                pA++;
            } 
            else if (pB < endB && (pA >= endA || b_ptr->col_indices[pB] < this->col_indices[pA])) {
                // 只有 B 有值
                C->col_indices.push_back(b_ptr->col_indices[pB]);

                ZZ tmp = b_ptr->values[pB] % q;
                if (tmp < 0) tmp += q;

                C->values.push_back(tmp);
                pB++;
            } 
            else {
                // 两者列索引相同，相加
                ZZ sum;
                add(sum, this->values[pA], b_ptr->values[pB]);
                if (sum != 0) { // 只有非零才存储
                    C->col_indices.push_back(this->col_indices[pA]);

                    ZZ tmp = sum % this->q;
                    if (tmp < 0) tmp += q;

                    C->values.push_back(tmp);
                }
                pA++; pB++;
            }
        }
        C->row_ptr[i+1] = C->col_indices.size();
    }
    C->sparsity = this->sparsity + B->getsparsity(); // 论文证明当中 sparsity 的扩张
    return C; 
}

std::unique_ptr<SparseMatrix> SparseMatrixCSR::operator*(const std::unique_ptr<SparseMatrix>& B) const {
    auto b_ptr = dynamic_cast<SparseMatrixCSR*>(B.get());
    if (!b_ptr) throw std::runtime_error("CSR required");

    long A_rows = this->rows;
    long B_cols = b_ptr->getcols();
    auto C = std::make_unique<SparseMatrixCSR>(A_rows, B_cols, this->q);
    
    C->row_ptr.resize(A_rows + 1);
    C->row_ptr[0] = 0;

    // 1. 我们需要为每一行独立存储计算出的非零元素，以避免线程冲突
    // 使用 vector 的 vector 来暂存每一行的结果
    std::vector<std::vector<long>> all_row_indices(A_rows);
    std::vector<std::vector<ZZ>> all_row_values(A_rows);

    // 2. 开启 OpenMP 并行计算行
    // schedule(dynamic) 适合稀疏矩阵，因为不同行的非零元数量差异很大
    #pragma omp parallel
    {
        // 每个线程私有的 SPA 结构，避免频繁申请内存
        std::vector<ZZ> spa_values(B_cols);
        std::vector<bool> spa_occupied(B_cols, false);
        std::vector<long> spa_indices;
        spa_indices.reserve(B_cols);

        #pragma omp for schedule(dynamic, 16)
        for (long i = 0; i < A_rows; ++i) {
            spa_indices.clear();

            // --- 核心计算逻辑 ---
            for (long j = this->row_ptr[i]; j < this->row_ptr[i+1]; ++j) {
                const ZZ& a_val = this->values[j];
                long a_col = this->col_indices[j];

                for (long k = b_ptr->row_ptr[a_col]; k < b_ptr->row_ptr[a_col+1]; ++k) {
                    long b_col = b_ptr->col_indices[k];
                    const ZZ& b_val = b_ptr->values[k];

                    if (!spa_occupied[b_col]) {
                        spa_occupied[b_col] = true;
                        spa_indices.push_back(b_col);
                        clear(spa_values[b_col]); 
                    }
                    
                    ZZ tmp;
                    mul(tmp, a_val, b_val);
                    add(spa_values[b_col], spa_values[b_col], tmp);
                }
            }

            std::sort(spa_indices.begin(), spa_indices.end());

            // 将结果存入当前行的私有暂存区
            for (long col : spa_indices) {
                if (this->q != 0) spa_values[col] %= this->q;

                if (spa_values[col] < 0) spa_values[col] += this->q; 
                
                if (spa_values[col] != 0) {
                    all_row_indices[i].push_back(col);
                    all_row_values[i].push_back(spa_values[col]);
                    clear(spa_values[col]); // 及时回收 NTL 内存
                }
                spa_occupied[col] = false;
            }
        }
    } // OpenMP 并行区结束

    // 3. 串行合并结果（将各行暂存区移动到 C 中）
    // 这一步虽然是串行，但主要是内存拷贝，速度很快
    for (long i = 0; i < A_rows; ++i) {
        C->col_indices.insert(C->col_indices.end(), 
                              all_row_indices[i].begin(), all_row_indices[i].end());
        C->values.insert(C->values.end(), 
                         std::make_move_iterator(all_row_values[i].begin()), 
                         std::make_move_iterator(all_row_values[i].end()));
        C->row_ptr[i+1] = C->col_indices.size();
    }

    C->sparsity = this->sparsity * B->getsparsity();
    return C;
}

std::unique_ptr<SparseMatrix> SparseMatrixCSR::addnewcolumn(const vec_ZZ& b) const {
    // 1. 维度检查
    if (b.length() != this->rows) {
        throw std::invalid_argument("Vector length must match matrix row count.");
    }

    // 2. 创建新的矩阵对象
    // 显式调用之前讨论过的构造函数初始化维度 (rows, cols + 1, k + 1, q)
    auto new_matrix = std::make_unique<SparseMatrixCSR>(this->rows, this->cols + 1, this->sparsity + 1, this->q);

    // 预估新容量以减少内存重分配
    // 新的非零元数量 = 原数量 + vector b 中的非零元数量
    long b_nnz = 0;
    for(long i = 0; i < b.length(); ++i) if (b[i] != 0) b_nnz++;
    
    new_matrix->values.reserve(this->values.size() + b_nnz);
    new_matrix->col_indices.reserve(this->col_indices.size() + b_nnz);
    new_matrix->row_ptr.resize(this->rows + 1);

    // 3. 执行合并逻辑
    new_matrix->row_ptr[0] = 0;
    long current_nnz = 0;

    for (long i = 0; i < this->rows; ++i) {
        // A. 拷贝原始矩阵第 i 行的所有元素
        for (long j = this->row_ptr[i]; j < this->row_ptr[i + 1]; ++j) {
            new_matrix->values.push_back(this->values[j]);
            new_matrix->col_indices.push_back(this->col_indices[j]);
            current_nnz++;
        }

        // B. 检查向量 b 的第 i 个元素
        // 如果不为 0，则作为该行最后一列（索引为原 cols）添加
        if (b[i] != 0) {
            new_matrix->values.push_back(b[i]);
            new_matrix->col_indices.push_back(this->cols); // 最后一列索引
            current_nnz++;
        }

        // C. 更新 row_ptr
        new_matrix->row_ptr[i + 1] = current_nnz;
    }

    // 返回 unique_ptr
    return new_matrix;
}

vec_ZZ generateSparseBernoulliVec(long n, const ZZ& q, double delta) {
    vec_ZZ v;
    v.SetLength(n);

    // 1. 计算概率 p = n^(-delta)
    double p = pow((double)n, -delta);

    // 2. 设置随机数引擎
    static random_device rd;
    static mt19937 gen(rd());
    
    // 伯努利分布：以 p 的概率返回 true
    bernoulli_distribution is_nonzero(p);
    long cnt = 0;
    for (long i = 0; i < n; ++i) {
        if (is_nonzero(gen)) {
            // 3. 如果不为 0，在 [1, q-1] 之间均匀采样
            // RandomBnd(q) 生成 [0, q-1]，我们需要排除 0
            cnt++;
            ZZ val;
            val = 0;
#ifdef WITH_ERROR
            do {
                val = RandomBnd(2);
            } while (val == 0); 
#endif
            v[i] = val;
        } else {
            // 4. 否则设为 0
            v[i] = 0;
        }
    }
    //cout << "Generating bernoulli vec with " << cnt << " non-zero items" << endl;
    return v;
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