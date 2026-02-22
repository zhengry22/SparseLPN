#pragma once
#include "polynomial.h"

Polynomial::Polynomial(std::vector<std::string> vars): variables(vars) {}
Polynomial::Polynomial() = default;

void Polynomial::add_term(double coeff, const Exponents& exps) {
        if (coeff == 0.0) return;
        if (exps.size() != variables.size()) {
            throw std::invalid_argument("Exponent vector size does not match number of variables");
        }
        terms[exps] += coeff;
}

// Get coefficient for a specific term
double Polynomial::get_coefficient(const Exponents& exps) const {
        auto it = terms.find(exps);
        if (it != terms.end()) return it->second;
        return 0.0;
    }

// Get total number of variables
int Polynomial::get_num_variables() const {
        return static_cast<int>(variables.size());
    }

// Get degree of polynomial (max sum of exponents)
int Polynomial::get_degree() const {
        int max_deg = 0;
        for (const auto& [exps, coeff] : terms) {
            int deg = 0;
            for (int e : exps) deg += e;
            max_deg = std::max(max_deg, deg);
        }
        return max_deg;
    }

// Polynomial addition
Polynomial Polynomial::operator+(const Polynomial& other) const {
        if (variables != other.variables) {
            throw std::invalid_argument("Polynomials have different variables");
        }
        Polynomial result(variables);
        result.terms = terms; // copy this polynomial
        for (const auto& [exps, coeff] : other.terms) {
            result.terms[exps] += coeff;
        }
        return result;
    }

// Polynomial multiplication
Polynomial Polynomial::operator*(const Polynomial& other) const {
        if (variables != other.variables) {
            throw std::invalid_argument("Polynomials have different variables");
        }
        Polynomial result(variables);
        for (const auto& [exps1, coeff1] : terms) {
            for (const auto& [exps2, coeff2] : other.terms) {
                Exponents new_exps(exps1.size());
                for (size_t i = 0; i < exps1.size(); ++i) {
                    new_exps[i] = exps1[i] + exps2[i];
                }
                result.terms[new_exps] += coeff1 * coeff2;
            }
        }
        return result;
    }

    // Print polynomial
void Polynomial::print() const {
    bool first = true;
    for (const auto& [exps, coeff] : terms) {
        if (!first && coeff >= 0) std::cout << "+";
        std::cout << coeff;
        for (size_t i = 0; i < variables.size(); ++i) {
            if (exps[i] != 0) std::cout << "*" << variables[i] << "^" << exps[i];
        }
        first = false;
    }
    if (first) std::cout << "0";
    std::cout << std::endl;
}