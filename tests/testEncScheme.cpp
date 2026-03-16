#include<iostream>
#include"EncryptionScheme.h"
#include"SparseMatrix.h"
#include"LinearHE.h"
#include <cstdlib>
using namespace std;
int main() {
    // Build up test here

    // 感觉参数的选取和多项式有关系

    long lambda = 2; 
    long tau = 1;
    auto lheptr = std::make_unique<Paillier>(lambda, tau);
    Polynomial poly;
    EncScheme myscheme(lambda, tau, std::move(lheptr), poly);
    
    // 建立了 EncScheme 之后
    KeyPair my_keypair = myscheme.keygen();
    myscheme.print();

    //ZZ mu1 = RandomBnd(myscheme.getmod());
    ZZ mu1 = RandomBnd(conv<ZZ>("10000000000"));
    cout << "plaintext 1: " << mu1 << endl;

    //ZZ mu2 = RandomBnd(myscheme.getmod());
    ZZ mu2 = RandomBnd(conv<ZZ>("10000000000"));
    cout << "plaintext 2: " << mu2 << endl;

    ZZ mu3 = RandomBnd(conv<ZZ>("10000000000"));
    cout << "plaintext 3: " << mu2 << endl;

    ZZ prod = mu1 * mu2;
    prod = prod % myscheme.getmod();
    prod = prod * mu3;
    prod = prod % myscheme.getmod();
    cout << "prod of plaintext: " << prod << endl;

    shescheme::Ciphertext ciphertext1 = myscheme.encrypt(my_keypair.secretkey, mu1);
    shescheme::Ciphertext expanded_ct1 =  myscheme.expand(my_keypair.secretkey, ciphertext1);

    shescheme::Ciphertext ciphertext2 = myscheme.encrypt(my_keypair.secretkey, mu2);
    shescheme::Ciphertext expanded_ct2 =  myscheme.expand(my_keypair.secretkey, ciphertext2);

    shescheme::Ciphertext ciphertext3 = myscheme.encrypt(my_keypair.secretkey, mu3);
    shescheme::Ciphertext expanded_ct3 =  myscheme.expand(my_keypair.secretkey, ciphertext3);

    auto expanded_ct = expanded_ct1 * expanded_ct2 * expanded_ct3;

    // 进行多项式运算

    ZZ compacted_ct = myscheme.compact(my_keypair.evalkey, my_keypair.secretkey, expanded_ct);
    ZZ decrypted_plaintext = myscheme.decrypt(compacted_ct);

    cout << "decrypted_plaintext: " << decrypted_plaintext << endl;
    ZZ difference = decrypted_plaintext - prod;
    double ratio = conv<double>(difference) / conv<double>(prod);
    cout << "The difference is: " << abs(difference) << endl;
    cout << "The difference ratio (diff / mu) is: " << abs(ratio) << endl;
    return 0;
}