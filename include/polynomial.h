#pragma once
#include<iostream>
#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <NTL/vec_ZZ.h>
#include <NTL/ZZ.h>
using namespace NTL;
using namespace std;
class Polynomial {
public:
    using Exponents = vector<int>;  // Each index corresponds to a variable
private:
    std::vector<std::string> variables;           // Variable names
    std::map<Exponents, ZZ> terms;            /*
                                                        a mapping vector, whose key is {exp1, exp2, ...} which
                                                        corresponds to the exponent of variables, and value is 
                                                        the coefficient of this term.

                                                        e.g. 
                                                            key --------- value
                                                        
                                                         {2, 4, 1} -------- 1.5

                                                         1.5(x_1)^2(x_2)^4(x_3)
                                                  */ 


public:
    // Constructor: define variables
    Polynomial(std::vector<std::string> vars);
    // Default constructor
    Polynomial();
    // Add a term: coefficient * x1^e1 * x2^e2 * ...
    void add_term(ZZ coeff, const Exponents& exps);
    // Get coefficient for a specific term
    ZZ get_coefficient(const Exponents& exps) const;
    // Get total number of variables
    int get_num_variables() const;
    // Get degree of polynomial (max sum of exponents)
    int get_degree() const;
    // Polynomial addition
    Polynomial operator+(const Polynomial& other) const;
    // Polynomial multiplication
    Polynomial operator*(const Polynomial& other) const;
    // Print polynomial
    void print() const;
};