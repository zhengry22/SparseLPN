#pragma once
#include "LinearHE.h"
#include <iostream>
#include <vector>
#include <Eigen/Sparse>
#include <NTL/ZZ_p.h>

#include <NTL/ZZ.h>
#include <iostream>
#include <algorithm>

using namespace std;
using namespace NTL;

class LHE {
// A Base class for everyone 
private:
    long lambda;
    long tau;
public:
    LHE(long lambda_, long tau_): lambda(lambda), tau(tau) {}
};


class Paillier: public LHE {
public:
    ZZ n, n_sq, g, lambda_priv, mu;

    /**
     * @brief 参数化密钥生成
     * @param lambda 安全参数 (如 128)
     * @param ell 明文长度要求 (比特数)
     */
    void keygen(long lambda, long ell) {
        long n_bit_len;

        // 1. 基于 lambda 映射到对应的 RSA/Paillier 模数强度 (NIST 标准)
        if (lambda >= 256)      n_bit_len = 15360;
        else if (lambda >= 192) n_bit_len = 7680;
        else if (lambda >= 128) n_bit_len = 3072;
        else if (lambda >= 112) n_bit_len = 2048;
        else                    n_bit_len = 1024; // 最小 1024

        // 2. 确保 n_bit_len 满足明文空间容量 ell
        // 因为 Paillier 的明文空间是 Z_n，所以 log2(n) 必须大于 ell
        if (n_bit_len < ell) {
            n_bit_len = ell + 1; 
        }

        cout << "[Paillier] Generating keys for lambda=" << lambda 
             << ", ell=" << ell << " (Total N bits: " << n_bit_len << ")" << endl;

        // 3. 生成大素数 p 和 q
        ZZ p, q;
        do {
            GenPrime(p, n_bit_len / 2);
            GenPrime(q, n_bit_len / 2);
        } while (GCD(p*q, (p-1)*(q-1)) != 1);
        
        n = p * q;
        n_sq = n * n;
        g = n + 1; // 优化：g = n + 1

        // 4. 计算私钥 lambda_priv = lcm(p-1, q-1)
        ZZ p_1 = p - 1;
        ZZ q_1 = q - 1;
        lambda_priv = (p_1 * q_1) / GCD(p_1, q_1); // 用这种方式计算最小公倍数

        // 5. 计算 mu = L(g^lambda_priv mod n^2)^-1 mod n
        ZZ u = PowerMod(g, lambda_priv, n_sq);
        ZZ l_u = (u - 1) / n;
        InvMod(mu, l_u, n);
    }

    // 加密
    ZZ encrypt(const ZZ& m) {
        ZZ r, c;
        do {
            RandomBnd(r, n);
        } while (GCD(r, n) != 1 || r == 0);

        // c = (1 + m*n) * r^n mod n^2
        ZZ gm = (1 + (m % n) * n) % n_sq;
        ZZ rn = PowerMod(r, n, n_sq);
        MulMod(c, gm, rn, n_sq);
        return c;
    }

    // 解密
    ZZ decrypt(const ZZ& c) {
        ZZ u = PowerMod(c, lambda_priv, n_sq);
        ZZ l_u = (u - 1) / n;
        ZZ m;
        MulMod(m, l_u, mu, n);
        return m;
    }

    // 同态加法
    ZZ add(const ZZ& c1, const ZZ& c2) {
        ZZ res;
        MulMod(res, c1, c2, n_sq);
        return res;
    }

    // 同态标量乘
    ZZ mul(const ZZ& c, const ZZ& k) {
        ZZ res;
        res = PowerMod(c, k, n_sq);
        return res;
    }
};