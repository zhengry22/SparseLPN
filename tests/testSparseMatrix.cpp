#pragma once
#include <iostream>
#include <vector>
#include <NTL/mat_ZZ.h>
#include <NTL/vec_ZZ.h>
#include "SparseMatrix.h"
#include <chrono>
//#define BASIC
//#define ALGEBRA
//#define MULT
#define PERFORMANCE
//#define NEW_COL
using namespace std;
using namespace NTL;

// 用于调试 SparseMatrix 是否能够正确的添加一列
    mat_ZZ removeLastColumn(const mat_ZZ& A) {
        long m = A.NumRows();
        long n = A.NumCols();
    
        if (n <= 0) {
            return mat_ZZ(); // 如果矩阵本来就没有列，返回空
        }
    
        mat_ZZ res;
        res.SetDims(m, n - 1); // 设置新维度为 m x (n-1)
    
        for (long i = 0; i < m; i++) {
            for (long j = 0; j < n - 1; j++) {
                res[i][j] = A[i][j];
            }
        }
    
        return res;
    }

int main() {
    // --- 1. 参数设置 ---
    long m = 16384; // 矩阵行数 (建议先用小规模验证逻辑)
    long n = 16384;  // 矩阵列数
    long p = 16384;
    long y = 170;
    int k1 = 40;   // 每行非零元
    int k2 = 40;
    int k3 = 30;
    ZZ q = conv<ZZ>("18446744073709551557");
    SparseMatrixCSRSampler sampler;
#ifdef BASIC
    // --- 2. 采样两个 RDiag 矩阵 ---
    // 假设采样器返回 std::unique_ptr
    
    auto uR1 = sampler.sample_RDiag(m, n, k1, q);
    auto uR2 = sampler.sample_RDiag(m, n, k2, q);

    // 为了方便操作，直接获取引用并向下转换为 CSR 类型
    auto& sA = dynamic_cast<SparseMatrixCSR&>(*uR1);
    auto& sB = dynamic_cast<SparseMatrixCSR&>(*uR2);

    cout << "--- 验证开始 ---" << endl;

    // --- 3. 验证 RDiag 结构定义 ---
    bool struct_ok = true;
    for (long i = 0; i < m; ++i) {
        long i_mod = (i + 1) % (n + 1); // 
        if (i_mod != 0) {
            long target_col = i_mod - 1; 
            bool found = false;
            // 遍历第 i 行的所有列索引
            for (long j = sA.row_ptr[i]; j < sA.row_ptr[i+1]; ++j) {
                if (sA.col_indices[j] == target_col) {
                    found = true; break;
                }
            }
            if (!found) {
                cout << "[结构错误] 第 " << i << " 行缺少对角元素，列: " << target_col << endl;
                struct_ok = false;
            }
        }
    }
    if (struct_ok) cout << "[PASS] RDiag 结构符合定义" << endl;

    // --- 4. 验证矩阵加法 ---
    // 手动包装指针以匹配你的接口：std::unique_ptr<SparseMatrix> operator+(const std::unique_ptr<SparseMatrix>&)
    auto ptr_for_add = std::make_unique<SparseMatrixCSR>(sB); 
    auto sSumBase = (*uR1) + uR2;
    auto& sSum = dynamic_cast<SparseMatrixCSR&>(*sSumBase);

    mat_ZZ dA = ToDense(sA);
    mat_ZZ dB = ToDense(sB);
    mat_ZZ dSum_gold = dA + dB;
    //cout << "dA: " << dA << endl;
    //cout << "dB: " << dB << endl;
    // 模运算
    for(long i=0; i<m; i++) for(long j=0; j<n; j++) dSum_gold[i][j] %= q;
    mat_ZZ todense = ToDense(sSum);
    //cout << "正确的矩阵加法: " << dSum_gold << endl;
    //cout << "实际做的加法: " << todense << endl;

    if (todense == dSum_gold) cout << "[PASS] 矩阵加法正确" << endl;
    else cout << "[FAIL] 矩阵加法不一致" << endl;

    // --- 5. 验证矩阵标量乘法 ---
    ZZ scalar = conv<ZZ>("10");
    auto sScalarBase = sA.operator*(scalar);
    auto& sScalar = dynamic_cast<SparseMatrixCSR&>(*sScalarBase);
    
    mat_ZZ dScalar_gold = dA * scalar;
    for(long i=0; i<m; i++) for(long j=0; j<n; j++) dScalar_gold[i][j] %= q;

    if (ToDense(sScalar) == dScalar_gold) cout << "[PASS] 标量乘法正确" << endl;
    else cout << "[FAIL] 标量乘法不一致" << endl;

    // --- 6. 验证矩阵-向量乘法 ---
    vec_ZZ v;
    v.SetLength(n);
    for(long i=0; i<n; i++) v[i] = RandomBnd(q);

    vec_ZZ sVecRes = sA.operator*(v);
    vec_ZZ dVecRes_gold = dA * v;
    for(long i=0; i<m; i++) dVecRes_gold[i] %= q;

    if (sVecRes == dVecRes_gold) cout << "[PASS] 矩阵向量乘法正确" << endl;
    else cout << "[FAIL] 矩阵向量乘法不一致" << endl;

    // --- 7. 验证矩阵-矩阵乘法 (SpGEMM) ---
    // 注意：sB_mul 的行数必须等于 sA 的列数
    auto uB_mul = sampler.sample_RDiag(n, p, k3, q);
    auto& sB_mul = dynamic_cast<SparseMatrixCSR&>(*uB_mul);
    
    auto ptr_for_mul = std::make_unique<SparseMatrixCSR>(sB_mul);
    auto sProdBase = (*uR1) * uB_mul;
    auto& sProd = dynamic_cast<SparseMatrixCSR&>(*sProdBase);
#ifdef MULT
    cout << "row number of result: " << sProd.getrows() << endl;
    cout << "col number of result: " << sProd.getcols() << endl;
#endif 
    mat_ZZ dB_mul = ToDense(sB_mul);
    mat_ZZ dProd_gold = dA * dB_mul;
#ifdef MULT
    cout << "A: " << dA << endl;
    cout << "B: " << dB_mul << endl;
#endif
    for(long i=0; i<m; i++) for(long j=0; j<p; j++) dProd_gold[i][j] %= q;
    mat_ZZ todenseprod = ToDense(sProd);
#ifdef MULT    
    cout << "正确的矩阵乘法: " << dProd_gold << endl;
    cout << "实际做的乘法: " << todenseprod << endl;
#endif
    if (todenseprod == dProd_gold) cout << "[PASS] 矩阵乘法(SpGEMM)正确" << endl;
    else cout << "[FAIL] 矩阵乘法(SpGEMM)不一致" << endl;
#endif
    cout << "接下来是大矩阵各种运算定律的正确性测试，" << endl;

#ifdef ALGEBRA    
    auto A = sampler.sample_RDiag(m, n, k1, q);
    auto B = sampler.sample_RDiag(m, n, k2, q);
    auto C = sampler.sample_RDiag(m, n, k3, q);
    auto D = sampler.sample_RDiag(n, p, k3, q);
    auto E = sampler.sample_RDiag(n, p, k2, q);
    auto F = sampler.sample_RDiag(p, y, k3, q);
    auto A_plus_A = (*A) + A;
    auto two_A = (*A) * conv<ZZ>("2");

    auto AA = dynamic_cast<SparseMatrixCSR&>(*A_plus_A);

    vec_ZZ x;
    x.SetLength(n);
    for(long i=0; i<n; i++) x[i] = RandomBnd(q);

    vec_ZZ res1 = (*A_plus_A) * x;
    vec_ZZ res2 = (*two_A) * x;

    // A + A = A * 2
    if ((*A_plus_A) == two_A) cout << "[PASS] A + A = A * 2 测试正确" << endl;
    else cout << "[FAIL] A + A = A * 2 测试错误" << endl;

    // (A + A) * x = (A * 2) * x
    if (res1 == res2) cout << "[PASS] (A + A) * x = (A * 2) * x 测试正确" << endl;
    else cout << "[FAIL] (A + A) * x = (A * 2) * x 测试错误" << endl;

    // (A + B) * x = A * x + B * x
    auto A_plus_B = (*A) + B;
    vec_ZZ res3 = (*A_plus_B) * x;
    vec_ZZ res4 = (*A) * x;
    vec_ZZ res5 = (*B) * x;
    vec_ZZ res6 = res4 + res5;
    for (long i = 0; i < res6.length(); ++i) {
        res6[i] %= q;
        // 如果需要保证结果为正数 (针对负数项)
        if (res6[i] < 0) res6[i] += q;
    }

    if (res3 == res6) cout << "[PASS] (A + B) * x = A * x + B * x 测试正确" << endl;
    else cout << "[FAIL] (A + B) * x = A * x + B * x 测试错误" << endl;

    // A + B = B + A
    auto B_plus_A = (*B) + A;
    if ((*A_plus_B) == B_plus_A) cout << "[PASS] A + B = B + A 测试正确" << endl;
    else cout << "[FAIL] A + B = B + A 测试错误" << endl;

    // (A + B) + C = A + (B + C)
    auto B_plus_C = (*B) + C;
    if (*((*A_plus_B) + C) == (*A) + (B_plus_C)) cout << "[PASS] A + (B + C) = A + (B + C) 测试正确" << endl;
    else cout << "[FAIL] A + (B + C) = A + (B + C) 测试错误" << endl;

    // (A * E) * F = A * (E * F)
    auto A_times_E = (*A) * E;
    auto E_times_F = (*E) * F;
    if (*((*A_times_E) * F) == (*A) * E_times_F) cout << "[PASS] (A * B) * C = A * (B * C) 测试正确" << endl;
    else cout << "[FAIL] (A * B) * C = A * (B * C) 测试错误" << endl;

    // A * (D + E) = A * D + A * E
    auto A_times_D = (*A) * D; 
    if (*((*A) * ((*D) + E)) == (*A_times_D) + A_times_E) cout << "[PASS] A * (D + E) = A * D + A * E 测试正确" << endl;
    else cout << "[FAIL] A * (D + E) = A * D + A * E 测试错误" << endl;

    // (A + B) * D = A * D + B * D
    auto B_times_D = (*B) * D;
    if (*((*A_plus_B) * D) == (*A_times_D) + B_times_D) cout << "[PASS]  A + (B + C) = A + (B + C) 测试正确" << endl;
    else cout << "[FAIL] A + (B + C) = A + (B + C) 测试错误" << endl; 
#endif    
#ifdef PERFORMANCE
    cout << "最后是性能测试。我们来比对 我们编写的系数矩阵乘法 和 一般的矩阵乘法 的 效率" << endl;
    auto M1 = sampler.sample_RDiag(m, n, k1, q);
    auto M2 = sampler.sample_RDiag(m, n, k2, q);
    
    auto start = std::chrono::high_resolution_clock::now();
    auto A1 = (*M1) * M2;
    auto end = std::chrono::high_resolution_clock::now();
    auto& AA = dynamic_cast<SparseMatrixCSR&>(*A1);
    std::chrono::duration<double> elapsed = end - start;
    cout << "第 1 轮乘法用时 " << elapsed.count() << " 秒, 含有 " << AA.values.size() << " 个元素" << endl;
    for (int i = 0; i < 16; i++) {
        auto M_tmp = sampler.sample_RDiag(m, n, k2, q);
        auto start = std::chrono::high_resolution_clock::now();
        A1 = (*A1) * M_tmp;
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        auto& AA = dynamic_cast<SparseMatrixCSR&>(*A1);
        cout << "第 " << i + 2 << " 轮乘法用时 " << elapsed.count() << " 秒, 含有 " << AA.values.size() << " 个元素" << endl;
    }
#endif   
#ifdef NEW_COL
    cout << "下面对于 插入新一列 进行测试: " << endl;
    auto A1 = sampler.sample_RDiag(n + 1, n, k1, q);
    auto AA1 = dynamic_cast<SparseMatrixCSR&>(*A1);
    auto b = generateSparseBernoulliVec(n + 1, q, 2);
    auto A1b = AA1.addnewcolumn(b);
    auto AA1b = dynamic_cast<SparseMatrixCSR&>(*A1b);

    auto A1dense = ToDense(AA1);
    auto A1b_dense = ToDense(AA1b);
    auto A1b_removecol = removeLastColumn(A1b_dense);

    if (A1dense == A1b_removecol) cout << "[PASS] 插入新一列测试成功" << endl;
    else cout << "[FAIL] 插入新一列测试失败" << endl;
#endif
    return 0;
}