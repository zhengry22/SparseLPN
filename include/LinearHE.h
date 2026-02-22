#include<iostream>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ.h>
using NTL::ZZ;
namespace lhescheme {

class EvaluationKey {
public:
    virtual ~EvaluationKey() = default;
};

class SecretKey {
public:
    virtual ~SecretKey() = default;
};

class PublicKey {
public:
    virtual ~PublicKey() = default;
};

}

class PaillierSecretKey: public lhescheme::SecretKey {
    ZZ lambda;
    ZZ mu;
public:
    PaillierSecretKey(ZZ lambda_, ZZ mu_);
};

class PaillierPublicKey: public lhescheme::PublicKey {
    ZZ n;
    ZZ g;
public:
    PaillierPublicKey(ZZ n_, ZZ g_);
};

class PaillierEvaluationKey: public lhescheme::EvaluationKey {
    // Its empty, since Paillier doesn't actually need Eval key
public:
    PaillierEvaluationKey();
};

class LHE {
// A Base class for everyone 
private:
    long lambda;
    long tau;
    lhescheme::EvaluationKey ek;
    lhescheme::SecretKey sk;
    lhescheme::PublicKey pk;
protected:
    void setEvalKey(lhescheme::EvaluationKey ek_);
    void setSecretKey(lhescheme::SecretKey sk_);
    void setPublicKey(lhescheme::PublicKey pk_);
public:
    LHE(long lambda, long tau);
};

class Paillier: public LHE {
public:
    ZZ n, n_sq, g, lambda_priv, mu;
    /**
     * @brief 参数化密钥生成
     * @param lambda 安全参数 (如 128)
     * @param ell 明文长度要求 (比特数)
     */
    using LHE::LHE;
    void keygen(long lambda, long ell);
    ZZ encrypt(const ZZ& m);
    ZZ decrypt(const ZZ& c);
    ZZ add(const ZZ& c1, const ZZ& c2);
    ZZ mul(const ZZ& c, const ZZ& k);
};
