#include <iostream>
#include <vector>
#include <NTL/mat_ZZ.h>
#include <NTL/vec_ZZ.h>
#include "SparseMatrix.h"
//#define MULT
using namespace std;
using namespace NTL;

int main() {
    // --- 1. 参数设置 ---
    long m = 100; // 矩阵行数 (建议先用小规模验证逻辑)
    long n = 200;  // 矩阵列数
    long p = 150;
    int k1 = 30;   // 每行非零元
    int k2 = 40;
    int k3 = 20;
    ZZ q = conv<ZZ>("65537");

    // --- 2. 采样两个 RDiag 矩阵 ---
    // 假设采样器返回 std::unique_ptr
    SparseMatrixCSRSampler sampler;
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

    return 0;
}