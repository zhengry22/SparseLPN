#include<iostream>
#include"LinearHE.h"
#include"polynomial.h"
#include<vector>
namespace shescheme {

// Structure of Ciphertext
struct Ciphertext {
    bool is_expanded; // 标记是普通向量还是矩阵形式
    std::vector<long long> data; // 向量或矩阵数据
    int depth; // 当前噪声深度
};

struct SecretKey {
    //TODO
};

class EvaluationKey {};

}

struct KeyPair {
    shescheme::SecretKey sk;
    lhescheme::EvaluationKey ek;
};

class EncScheme {
    /*
        Evaluate one polynomial per scheme; We should obtain the polynomial first,
        and then the parameters
    */ 
private:
    // The 3 parameters below are independent of other parameters

    // The polynomial to evaluate (It defines the security parameters) 
    Polynomial poly;

    // Security Standards
    int lambda;    // Safety Parameter 
    int tau;       // Capacity 

    // LPN parameters, which is initialized using data of poly, lambda and tau.
    double delta;  // Noise rate
    long long q;   // Modulus
    long long n;   // Number of Columns of A, which could be calculated using lamdba and tau
    int k;         // Sparsity
    
    LHE lhe;       // Linear HE

    int get_poly_degree() const;
    static double generate_valid_delta(int lambda_, int tau_);
    static long long generate_valid_q(int lambda_, int tau_, const Polynomial& poly);
    static long long generate_valid_n(int lambda_, int tau_, int delta); // Need delta
    static int generate_valid_k(int lambda_, int tau_, int n);

public:
    /**
     * @param lambda_ Safety parameter
     * @param max_depth How many monomials
     */
    EncScheme(int lambda_, int items, LHE lhe_, Polynomial poly_);

    // --- Core Operations ---


    KeyPair keygen(const shescheme::Ciphertext& ct, const lhescheme::EvaluationKey& ek);
    shescheme::Ciphertext encrypt(long long mu);
    shescheme::Ciphertext expand(const shescheme::Ciphertext& ct);
    shescheme::Ciphertext add(const shescheme::Ciphertext& ct1, const shescheme::Ciphertext& ct2);
    shescheme::Ciphertext mult(const shescheme::Ciphertext& ct1, const shescheme::Ciphertext& ct2);
    shescheme::Ciphertext compact(const shescheme::Ciphertext& expanded_ct, const shescheme::Ciphertext& ek);
    long long decrypt(const shescheme::Ciphertext& ct);
};

class BatchEncScheme: EncScheme {
    // To evaluate a batch of polynomials and various variants
};