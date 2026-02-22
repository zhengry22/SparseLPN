#include <iostream>
#include <chrono>
#include <cassert>
#include <NTL/ZZ.h>
#include "LinearHE.h" // 假设你的 Paillier 类定义在此头文件中

using namespace std;
using namespace NTL;
using NTL::ZZ;
#define TEST_PAILLIER

// 简单的辅助函数，用于打印测试结果
void print_test_result(const string& test_name, bool success) {
    cout << (success ? "[  OK  ] " : "[FAILED] ") << test_name << endl;
}

int main() {
#ifdef TEST_PAILLIER
cout << "--- Starting Advanced Paillier LHE Unit Test ---" << endl;

    // 1. 系统初始化
    long lambda = 128; 
    long ell = 100;    
    Paillier lhe(lambda, 0);

    auto start = chrono::high_resolution_clock::now();
    lhe.keygen(lambda, ell);
    auto end = chrono::high_resolution_clock::now();
    
    chrono::duration<double> diff = end - start;
    cout << "Keygen Time: " << diff.count() << "s | N size: " << NumBits(lhe.n) << " bits" << endl;

    // ---------------------------------------------------------
    // 2. 随机压力测试 (100轮)
    // ---------------------------------------------------------
    bool all_passed = true;
    for(int i = 0; i < 100; ++i) {
        ZZ m = RandomBnd(lhe.n); // 生成 [0, n-1] 之间的随机数
        ZZ ct = lhe.encrypt(m);
        if (lhe.decrypt(ct) != m) {
            all_passed = false;
            break;
        }
    }
    print_test_result("Random Stress Test (100 rounds)", all_passed);

    // ---------------------------------------------------------
    // 3. 边界值测试 (最大值与最小值)
    // ---------------------------------------------------------
    ZZ m_max = lhe.n - 1;
    ZZ m_min = conv<ZZ>("1");
    bool boundary_ok = (lhe.decrypt(lhe.encrypt(m_max)) == m_max) && 
                       (lhe.decrypt(lhe.encrypt(m_min)) == m_min);
    print_test_result("Boundary Values (1 and n-1)", boundary_ok);

    // ---------------------------------------------------------
    // 4. 同态累加测试 (Σ m_i)
    // ---------------------------------------------------------
    ZZ expected_sum = conv<ZZ>("0");
    ZZ combined_ct = lhe.encrypt(expected_sum);
    for(int i = 1; i <= 10; ++i) {
        ZZ m_i = conv<ZZ>(i * 100);
        expected_sum += m_i;
        ZZ ct_i = lhe.encrypt(m_i);
        combined_ct = lhe.add(combined_ct, ct_i);
    }
    print_test_result("Homomorphic Summation (10 items)", lhe.decrypt(combined_ct) == expected_sum);

    // ---------------------------------------------------------
    // 5. 标量乘法与加法结合测试: a*m1 + b*m2
    // ---------------------------------------------------------
    ZZ ma = conv<ZZ>("50"), mb = conv<ZZ>("100");
    ZZ ka = conv<ZZ>("5"), kb = conv<ZZ>("10");
    // 计算 E(50)*5 + E(100)*10 = E(50*5 + 100*10) = E(1250)
    ZZ ct_a = lhe.encrypt(ma);
    ZZ ct_b = lhe.encrypt(mb);
    ZZ res_ct = lhe.add(lhe.mul(ct_a, ka), lhe.mul(ct_b, kb));
    
    ZZ expected_combined = (ma * ka) + (mb * kb);
    print_test_result("Linear Combination (a*m1 + b*m2)", lhe.decrypt(res_ct) == expected_combined);

    // ---------------------------------------------------------
    // 6. 溢出行为测试 (解密结果应符合模 n 运算)
    // ---------------------------------------------------------
    // 在 Paillier 中，明文空间是 Z_n。如果结果 > n，则会发生回绕
    ZZ m_large = lhe.n - 10;
    ZZ k_overflow = conv<ZZ>("20");
    ZZ overflow_ct = lhe.mul(lhe.encrypt(m_large), k_overflow);
    ZZ expected_overflow = (m_large * k_overflow) % lhe.n;
    
    print_test_result("Plaintext Space Overflow (Result mod n)", lhe.decrypt(overflow_ct) == expected_overflow);

    cout << "--- All Detailed Tests Completed ---" << endl;
#endif
    return 0;
}