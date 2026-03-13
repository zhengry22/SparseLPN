#include"EncryptionScheme.h"
#include <math.h>
#include <random>
#include <cmath>
//#define TEST_EXPAND
//#define TEST_COMPACT
//#define TEST_CORRECT
namespace shescheme {

}

void print_NTL_vec(const vec_ZZ& vec) {
    for (auto &e:vec) {
        cout << e << " ";
    }
    cout << endl;
}

template<typename T>
void print_vec(const vector<T>& vec) {
    for (auto &e:vec) {
        cout << e << " ";
    }
    cout << endl;
}

void print_dense_matrix(const mat_ZZ& mat) {
    for (long i = 0; i < mat.NumRows(); i++) {
        for (long j = 0; j < mat.NumCols(); j++) {
            cout << mat[i][j] << " ";
        }
        cout << endl;
    }
}

int EncScheme::get_poly_degree() const {
    return this->poly.get_degree();  
}

double EncScheme::generate_valid_delta(long lambda_, long tau_){
    /*
        Use a number in (0, 1) to initialize
    */    
    constexpr double delta = 0.8;
    if ((delta >= 1) || (delta) <= 0) {throw std::invalid_argument("delta should be in (0, 1)");}
    cout << "generate delta: " << delta << endl; 
    return delta;
}

ZZ EncScheme::generate_valid_q(
    long lambda_,
    long tau_,
    const Polynomial& poly
) {
    //constexpr double alpha = 0.1;
    //constexpr double beta  = 0.8;
    // TODO
    // 先随机生成一个比较小的 mod
    //cout << "generate q: " << conv<ZZ>("1073741827") << endl;
#ifdef TEST_CORRECT
    return conv<ZZ>("17");
#endif
    return conv<ZZ>("65537");
}

// Private functions
long long EncScheme::generate_valid_n(long lambda_, long tau_, double delta){
    constexpr long long c = 7;
    if (tau_ == 1) {
        //cout << "generate n: " << 6712 << endl;
        return 2;
    }
    auto ret = static_cast<long long>(
        std::ceil(std::pow(tau_, (2 / delta)) * std::pow(lambda_, (c / delta))) // Use ceiling to make sure it does not get smaller
    );
    cout << "generate n: " << ret << endl;
    return ret;
}

long EncScheme::generate_valid_k(long lambda_, long tau_, long long n){
    /*
        Generate a sparsity-index k that is:
            1. (much) smaller than n
            2. k <= \sqrt{log{\tau} - 1}
            3. is an odd number  
    */ 
    // Handle special cases here
    if (tau_ == 1) {
        //cout << "generate k: " << 2 << endl;
        return 1; // 我们可以先考虑使用平凡的 k = 2 作为 sparsity 看一下效果
    } 

    bool my_random = false;

    double log_tau = std::log(static_cast<double>(tau_));
    int k_max = static_cast<int>(std::floor(std::sqrt(log_tau - 1)));

    // k_max is at most n / 2
    k_max = std::min(k_max, static_cast<int>(n / 2));

    // k_max is at least 1
    k_max = std::max(k_max, 1);

    // k_max is odd
    if (k_max % 2 == 0) k_max -= 1;
    if (k_max < 1) k_max = 1;

    int k = k_max;
    if (my_random) { // Generate a random odd number
        int num_odd = (k_max + 1) / 2; // # of odd numbers
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(0, num_odd - 1);
        int r = dis(gen);
        k = 1 + 2 * r; // odd number between [1, k_max]
    }

    cout << "generate k: " << k << endl;
    return k;
}

EncScheme::EncScheme(long lambda_, long items, std::unique_ptr<LHE> lhe_, Polynomial poly_): 
    lambda(lambda_), tau(items), lhe(std::move(lhe_)), poly(poly_),
    delta(generate_valid_delta(lambda_, items)), 
    q(generate_valid_q(lambda_, items, poly_)), 
    n(generate_valid_n(lambda_, items, delta)), 
    k(generate_valid_k(lambda_, items, n)){    

    // Check whether the parameters are valid:
    if (k >= n) throw std::invalid_argument("k >= n!");
}

void EncScheme::print() {
    // Print out all the parameters
    std::cout << "lambda: " << this->lambda << std::endl;
    std::cout << "tau: " << this->tau << std::endl;
    std::cout << "delta: " << this->delta << std::endl;
    std::cout << "q: " << this->q << std::endl;
    std::cout << "n: " << this->n << std::endl;
    std::cout << "k: " << this->k << std::endl;
}

//-------------------Below are the functions used for the SHE scheme----------------------

vec_ZZ sampleZZVector(long n, const ZZ& q) {
    vec_ZZ v;
    v.SetLength(n);

    for (long i = 0; i < n; i++) {
        // RandomBnd(q) 返回一个 [0, q-1] 之间的随机大整数
        v[i] = RandomBnd(q);
    }
    return v;
}

std::unique_ptr<SparseMatrix> EncScheme::GSWEnc(const vec_ZZ& s, ZZ& mu) {
    if (s.length() != this->n) {
        throw std::invalid_argument("[EncScheme::GSWEnc] s must have a length of n! ");
    }
    long long l = this->n + 1;
    
    auto AA = this->sampler.sample_RDiag(l, n, k, q);
    auto& A = dynamic_cast<SparseMatrixCSR&>(*AA);

    vec_ZZ e = generateSparseBernoulliVec(l, this->q, this->delta);
    vec_ZZ s_tilde = -s;
    s_tilde.append(conv<ZZ>("1"));

    auto b = A * s + e + mu * s_tilde;
    if (b.length() != this->n + 1) {
        throw std::invalid_argument("[EncScheme::GSWEnc] b is not n + 1 long! ");
    }
    for(auto &e: b) {
        e %= this->q;
        if (e < 0) e += q;
    } 
    auto ret = A.addnewcolumn(b);
    if (ret->getrows() != A.getrows()) {
        throw std::invalid_argument("[EncScheme::GSWEnc] ret's rows doesn't = A's rows! ");
    }
    if (ret->getcols() != A.getcols() + 1) {
        throw std::invalid_argument("[EncScheme::GSWEnc] ret's cols doesn't = A's cols + 1! ");
    }
    return ret;
}

KeyPair EncScheme::keygen() {
    cout << "calling keygen..." << endl;
    // Note that in Paillier, there is no Eval Key, and we also need the public key
    ZZ newmod = (*lhe).keygen(this->lambda, this->tau);

    // If we are using Paillier, then the modulus of the scheme must be changed
    this->q = newmod;

    // Generate keys
    auto raw_evalkey = (*lhe).getEvalKey();
    auto raw_publickey = (*lhe).getPublicKey();
    auto raw_secretkey = (*lhe).getSecretKey();

    // Generate s and t
    vec_ZZ s, t;
    s = sampleZZVector(this->n, this->q);
    t = sampleZZVector(this->n, this->q);
#ifdef TEST_CORRECT
    cout << "generate s: ";
    print_NTL_vec(s);
    cout << "generate t: ";
    print_NTL_vec(t);
#endif 
    // vec_ZZ s_tilde = -s;
    // vec_ZZ t_tilde = -t;
    // s_tilde.append(conv<ZZ>("1"));
    // t_tilde.append(conv<ZZ>("1"));

    // vector<std::unique_ptr<SparseMatrix>> C_ek;
    // vec_ZZ vec_ct;
    // for (long long i = 0; i < (n + 1); i++) {
    //     C_ek.push_back(GSWEnc(s, t_tilde[i]));
    //     vec_ct.append(this->lhe->encrypt(s_tilde[i]));
    // }

    shescheme::EvaluationKey ev_key;
    if (raw_evalkey) {
        ev_key.ek_lhe = std::make_unique<lhescheme::EvaluationKey>(*raw_evalkey);
    }
    // ev_key.C_ek = std::move(C_ek);
    // ev_key.vec_ct = vec_ct;
 
    shescheme::SecretKey sc_key;
    if (raw_secretkey) {
        sc_key.sk_lhe = std::make_unique<lhescheme::SecretKey>(*raw_secretkey);
    }
    sc_key.s = s;
    sc_key.t = t;

    // Finally generate keypair and return
    KeyPair keypair = { std::move(ev_key), std::move(sc_key) };
    return keypair;
}

shescheme::Ciphertext EncScheme::encrypt(const shescheme::SecretKey &sk, ZZ& mu) {
    cout << "encrypting... " << endl;
    auto aa = this->sampler.sample_RDiag(1, this->n, this->k, this->q);
    auto& a = dynamic_cast<SparseMatrixCSR&>(*aa);
#ifdef TEST_CORRECT
    cout << "generate a: ";
    print_dense_matrix(ToDense(a));
#endif
    if (a.getrows() != 1) {
        throw std::invalid_argument("[Encscheme::encrypt] a must have only 1 row! ");
    }
    if (a.row_ptr.size() != 2) {
        throw std::invalid_argument("[Encscheme::encrypt] invalid a row_ptr size! ");
    }
    if (a.getcols() != sk.t.length()) {
        throw std::invalid_argument("[Encscheme::encrypt] a and t not the same size! ");
    }

    vec_ZZ e = generateSparseBernoulliVec(1, this->q, this->delta);
    if (e.length() != 1) {
        throw std::invalid_argument("[Encscheme::encrypt] e is not a scalar! ");
    }
#ifdef TEST_CORRECT
    cout << "generate e: "; 
    print_NTL_vec(e);
#endif 
    vec_ZZ m; 
    m.append(mu);
#ifdef TEST_CORRECT
    cout << "t vector: " ;
    print_NTL_vec(sk.t);
    cout << "generate m: ";
    print_NTL_vec(m);
#endif 
    auto prod = (a * sk.t) + e + m;
    auto c = a.addnewcolumn(prod);
    auto& cc = dynamic_cast<SparseMatrixCSR&>(*c);
#ifdef TEST_CORRECT
    cout << "product: ";
    print_NTL_vec(prod);
    auto dense = ToDense(cc);
    cout << "print c: ";
    print_dense_matrix(dense);
#endif
    if (c->getrows() != 1) {
        throw std::invalid_argument("[Encscheme::encrypt] c must have only 1 row! ");
    }
    if (cc.row_ptr.size() != 2) {
        throw std::invalid_argument("[Encscheme::encrypt] invalid c row_ptr size! ");
    }

    shescheme::Ciphertext ct = {false, std::move(c), 1};
    return ct;
}

shescheme::Ciphertext EncScheme::expand(const shescheme::SecretKey &sk, const shescheme::Ciphertext& ct) {
    cout << "calling expand... " << endl;
    if (ct.data->getrows() != 1) {
        throw std::invalid_argument("[Encscheme::expand] ct must have only 1 row! ");
    }
    if (ct.is_expanded) {
        throw std::invalid_argument("[Encscheme::expand] ct must be unexpanded! ");
    }

    // 为了调试方便，目前暂且只考虑 CSR 的情况，我们动态转换一下看看是不是 CSR 格式的稀疏矩阵
    auto* csr_ptr = dynamic_cast<SparseMatrixCSR*>(ct.data.get());

    if (csr_ptr) {
        // 转换成功：它确实是一个 SparseMatrixCSR 实例
        std::cout << "This is a SparseMatrixCSR." << std::endl;
        if (csr_ptr->row_ptr.size() != 2) {
            throw std::invalid_argument("[Encscheme::expand] invalid ct row_ptr size! ");
        }
        // 你现在可以访问 csr_ptr-> 特有的 CSR 成员了
    } else {
        // 转换失败：它是基类或其他派生类
        std::cout << "This is NOT a SparseMatrixCSR." << std::endl;
    }
    
    vec_ZZ c_vec = ct.data->getRowAsVec(-1);
    if (c_vec.length() != this->n + 1) {
        cout << "c_vec.length: " << c_vec.length() << endl;
        throw std::invalid_argument("[Encscheme::expand] c_vec does not have length n + 1! ");
    }
#ifdef TEST_CORRECT
    cout << "ciphertext c: ";
    print_NTL_vec(c_vec);
#endif

    // 为了内存优化，改成在这里计算 C_{ek_i}. 首先我们需要计算 s_tilde, t_tilde
    vec_ZZ s_tilde = -sk.s;
    vec_ZZ t_tilde = -sk.t;
    s_tilde.append(conv<ZZ>("1"));
    t_tilde.append(conv<ZZ>("1"));
#ifdef TEST_CORRECT
    cout << "s_tilde: ";
    print_NTL_vec(s_tilde);
    cout << "t_tilde: ";
    print_NTL_vec(t_tilde);
#endif     
    cout << "Adding matrices together in expand! " << endl;
    
    auto M = GSWEnc(sk.s, t_tilde[0]);
    auto C = *M * c_vec[0];

#ifdef TEST_CORRECT
    auto& M0 = dynamic_cast<SparseMatrixCSR&>(*M);
    cout << "M0: " << endl;
    print_dense_matrix(ToDense(M0));
    auto& C0 = dynamic_cast<SparseMatrixCSR&>(*C);
    cout << "C0: " << endl;
    print_dense_matrix(ToDense(C0));
#endif

    for (long i = 1; i < this->n + 1; i++) {
        auto MM = GSWEnc(sk.s, t_tilde[i]);
        C = (*C) + (*MM * c_vec[i]);

#ifdef TEST_CORRECT
        auto& Mi = dynamic_cast<SparseMatrixCSR&>(*MM);
        cout << "M" << i << ": " << endl;
        print_dense_matrix(ToDense(Mi));
        auto& Ci = dynamic_cast<SparseMatrixCSR&>(*C);
        cout << "C" << i << ": " << endl;
        print_dense_matrix(ToDense(Ci));
#endif

#ifdef TEST_EXPAND
        cout << "Finish adding: i = " << i << endl; 
#endif  
    }
    //cout << "[Encscheme::expand] finished adding" << endl;
    shescheme::Ciphertext ciphertext = {true, std::move(C), 1};
    //cout << "[Encscheme::expand] Returning Ciphertext..." << endl;
    return ciphertext;
}

ZZ EncScheme::compact(const shescheme::EvaluationKey& ek, const shescheme::SecretKey &sk, const shescheme::Ciphertext& ct) {
    cout << "[Encscheme::compact] calling compact!" << endl;
    if (!ct.is_expanded) {
        throw std::invalid_argument("[Encscheme::compact] ct is not expanded");
    }

    vec_ZZ c = ct.data->getRowAsVec(-1); // last row of C

    if (c.length() != this->n + 1) {
        throw std::invalid_argument("[Encscheme::compact] lenght of c is not n + 1! ");
    }

    // LHE eval
    auto lhe_ek = ek.ek_lhe.get();
    // 为了内存优化，改成在这里计算 C_{ek_i}. 首先我们需要计算 s_tilde, t_tilde
    vec_ZZ s_tilde = -sk.s;
    vec_ZZ t_tilde = -sk.t;
    s_tilde.append(conv<ZZ>("1"));
    t_tilde.append(conv<ZZ>("1"));
    //cout << "[Encscheme::compact] using lhe in compact..." << endl;
    ZZ ct0 = lhe->encrypt(s_tilde[0]);
    //cout << "[Encscheme::compact] ready to enter loop..." << endl;
    ZZ ret = lhe->mul(ct0, c[0], *lhe_ek);
    for (long long i = 1; i < this->n + 1; i++) {
        ZZ ct = lhe->encrypt(s_tilde[i]);
        ret = lhe->add(ret, lhe->mul(ct, c[i], *lhe_ek), *lhe_ek); 
#ifdef TEST_COMPACT
        cout << "Finish adding: i = " << i << endl; 
#endif
    }
    //cout << "[Encscheme::compact] ready to return..." << endl;
    return ret;
}

ZZ EncScheme::decrypt(const ZZ& ct) {
    // Simply decrypt
    //cout << "[Encscheme::decrypt] calling decrypt! " << endl;
    ZZ mu = lhe->decrypt(ct);
    //cout << "[Encscheme::decrypt] returning in decrypt" << endl;
    ZZ ret = mu % this->q;
    if (ret < 0) ret += q;
    return ret;
}