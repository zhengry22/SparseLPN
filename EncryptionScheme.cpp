#pragma once
#include"EncryptionScheme.h"
#include <math.h>
#include <random>

namespace shescheme {

}

int EncScheme::get_poly_degree() const {
    return this->poly.get_degree();  
}

double EncScheme::generate_valid_delta(int lambda_, int tau_){
    /*
        Use a number in (0, 1) to initialize
    */    
    constexpr double delta = 0.4;
    if ((delta >= 1) || (delta) <= 0) {throw std::invalid_argument("delta should be in (0, 1)");}
    return delta;
}

long long EncScheme::generate_valid_q(
    int lambda_,
    int tau_,
    const Polynomial& poly
) {
    constexpr double alpha = 0.1;
    constexpr double beta  = 0.8;

    int d = poly.get_degree();  

    return static_cast<long long>(
        std::ceil(lambda_ + alpha * std::sqrt(tau_) + beta * d) // Use ceiling to make sure it does not get smaller
    );
}

// Private functions
long long EncScheme::generate_valid_n(int lambda_, int tau_, int delta){
    //TODO

    constexpr long long c = 10;
    return static_cast<long long>(
        std::ceil(std::pow(tau_, (2 / delta)) * std::pow(lambda_, (c / delta))) // Use ceiling to make sure it does not get smaller
    );
}

int EncScheme::generate_valid_k(int lambda_, int tau_, int n){
    /*
        Generate a sparsity-index k that is:
            1. (much) smaller than n
            2. k <= \sqrt{log{\tau} - 1}
            3. is an odd number  
    */ 
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

    return k;
}

EncScheme::EncScheme(int lambda_, int items, LHE lhe_, Polynomial poly_): 
    lambda(lambda_), tau(items), lhe(lhe_), poly(poly_),
    delta(generate_valid_delta(lambda_, items)), 
    q(generate_valid_q(lambda_, items, poly_)), 
    n(generate_valid_n(lambda_, items, delta)), 
    k(generate_valid_k(lambda_, items, n)){    

    // Check whether the parameters are valid:
    if (k >= n) throw std::invalid_argument("k >= n!");
}

//-------------------Below are the functions used for the SHE scheme----------------------

