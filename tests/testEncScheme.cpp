#include<iostream>
#include"EncryptionScheme.h"
#include"SparseMatrix.h"
#include"LinearHE.h"
using namespace std;

int main() {
    // Build up test here
    long lambda = 2; 
    long tau = 5;
    auto lheptr = std::make_unique<Paillier>(lambda, tau);
    Polynomial poly;
    EncScheme myscheme(lambda, tau, std::move(lheptr), poly);
    // 建立了 EncScheme 之后
    KeyPair my_keypair = myscheme.keygen();
    ZZ mu = RandomBnd(myscheme.getmod());
    cout << "plaintext: " << mu << endl;

    shescheme::Ciphertext ciphertext = myscheme.encrypt(my_keypair.secretkey, mu);
    shescheme::Ciphertext expanded_ct =  myscheme.expand(my_keypair.evalkey, ciphertext);
    ZZ compacted_ct = myscheme.compact(my_keypair.evalkey, expanded_ct);
    ZZ decrypted_plaintext = myscheme.decrypt(compacted_ct);

    cout << "decrypted_plaintext" << decrypted_plaintext << endl;
    return 0;
}