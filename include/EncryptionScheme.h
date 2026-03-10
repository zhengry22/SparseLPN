#pragma once
#include<iostream>
#include"LinearHE.h"
#include"SparseMatrix.h"
#include"polynomial.h"
#include<vector>
#include <algorithm>
#include<memory>
namespace shescheme {

// Structure of Ciphertext
struct Ciphertext {
public:
    bool is_expanded; // 标记是普通向量还是矩阵形式
    std::unique_ptr<SparseMatrix> data; // 向量或矩阵数据
    int depth; // 当前噪声深度

    shescheme::Ciphertext operator+(const shescheme::Ciphertext& ct) {
        if (this->is_expanded != ct.is_expanded) {
            throw std::invalid_argument("[Ciphertext::Add] can't add non-expand with expand! ");
        }
        if (this->data->getrows() != ct.data->getrows()) {
            throw std::invalid_argument("[Ciphertext::Add] row number doesn't match! ");
        }
        if (this->data->getcols() != ct.data->getcols()) {
            throw std::invalid_argument("[Ciphertext::Add] column number doesn't match! ");
        }
        if (this->data->getmod() != ct.data->getmod()) {
            throw std::invalid_argument("[Ciphertext::Add] mod q doesn't match! ");
        }
        shescheme::Ciphertext ret;
        ret.is_expanded = this->is_expanded;
        ret.data = *(this->data) + ct.data;
        // TODO 目前还不太明白乘法深度是怎么弄的
        ret.depth = max(this->depth, ct.depth);
        return ret;
    }

    shescheme::Ciphertext operator*(const shescheme::Ciphertext& ct) {
        if (this->is_expanded != ct.is_expanded) {
            throw std::invalid_argument("[Ciphertext::Mult] can't mult non-expand with expand! ");
        }
        if (this->data->getcols() != ct.data->getrows()) {
            throw std::invalid_argument("[Ciphertext::Mult] can't mult these 2 matrices! ");
        }
        if (this->data->getmod() != ct.data->getmod()) {
            throw std::invalid_argument("[Ciphertext::Mult] mod q doesn't match! ");
        }
        shescheme::Ciphertext ret;
        ret.is_expanded = this->is_expanded;
        ret.data = *(this->data) * ct.data;
        // TODO 目前还不太明白乘法深度是怎么弄的
        ret.depth = this->depth + ct.depth;
        return ret;
    }
};

struct SecretKey {
    // 含有 s, t, sk_LHE.
    vec_ZZ s;
    vec_ZZ t;
    std::unique_ptr<lhescheme::SecretKey> sk_lhe;
};

struct EvaluationKey {
    /**
     * @param ek_lhe eval key
     * @param C_ek matrices
     * @param vec_ct vector
     */
    // 含有 ek_LHE, C_ek1, ... C_ek{n + 1}, ct1, ..., ct_l.
    std::unique_ptr<lhescheme::EvaluationKey> ek_lhe;
    vector<std::unique_ptr<SparseMatrix>> C_ek;
    vec_ZZ vec_ct;
};

}

struct KeyPair {
    shescheme::EvaluationKey evalkey;
    shescheme::SecretKey secretkey;
};

class EncScheme {
    /*
        Evaluate one polynomial per scheme; We should obtain the polynomial first,
        and then the parameters
    */ 
private:
    // The 3 parameters below are independent of other parameters

    // The polynomial to evaluate (It defines the security parameters) 
    
    SparseMatrixCSRSampler sampler;
    // Security Standards
    long lambda;    // Safety Parameter 
    long tau;       // Capacity 
    std::unique_ptr<LHE> lhe;       // Linear HE
    Polynomial poly;
    // LPN parameters, which is initialized using data of poly, lambda and tau.
    double delta;  // Noise rate
    ZZ q;   // Modulus
    long long n;   // Number of Columns of A, which could be calculated using lamdba and tau
    long k;         // Sparsity

    int get_poly_degree() const;
    static double generate_valid_delta(long lambda_, long tau_);
    static ZZ generate_valid_q(long lambda_, long tau_, const Polynomial& poly);
    static long long generate_valid_n(long lambda_, long tau_, double delta); // Need delta
    static long generate_valid_k(long lambda_, long tau_, long long n);

public:
    /**
     * @param lambda_ Safety parameter
     * @param max_depth How many monomials
     */
    EncScheme(long lambda_, long items, std::unique_ptr<LHE> lhe_, Polynomial poly_);

    // --- Core Operations ---

    std::unique_ptr<SparseMatrix> GSWEnc(const vec_ZZ& s, ZZ& mu);
    KeyPair keygen();
    shescheme::Ciphertext encrypt(const shescheme::SecretKey &sk, ZZ& mu);
    // shescheme::Ciphertext expand(const shescheme::EvaluationKey &ek, const shescheme::Ciphertext& ct);
    // ZZ compact(const shescheme::EvaluationKey& ek, const shescheme::Ciphertext& ct);
    shescheme::Ciphertext expand(const shescheme::SecretKey &sk, const shescheme::Ciphertext& ct);
    ZZ compact(const shescheme::EvaluationKey& ek, const shescheme::SecretKey &sk, const shescheme::Ciphertext& ct);
    ZZ decrypt(const ZZ& ct);
    ZZ getmod() {return this->q;}
    // For debugging
    void print();
};

class BatchEncScheme: EncScheme {
    // To evaluate a batch of polynomials and various variants
};