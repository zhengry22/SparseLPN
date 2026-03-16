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

    ZZ mu = RandomBnd(myscheme.getmod());
    cout << "plaintext: " << mu << endl;

    shescheme::Ciphertext ciphertext = myscheme.encrypt(my_keypair.secretkey, mu);
    shescheme::Ciphertext expanded_ct =  myscheme.expand(my_keypair.secretkey, ciphertext);

    // 进行多项式运算

    ZZ compacted_ct = myscheme.compact(my_keypair.evalkey, my_keypair.secretkey, expanded_ct);
    ZZ decrypted_plaintext = myscheme.decrypt(compacted_ct);

    cout << "decrypted_plaintext: " << decrypted_plaintext << endl;
    ZZ difference = decrypted_plaintext - mu;
    double ratio = conv<double>(difference) / conv<double>(mu);
    cout << "The difference is: " << abs(difference) << endl;
    cout << "The difference ratio (diff / mu) is: " << abs(ratio) << endl;
    return 0;
}