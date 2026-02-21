#include<iostream>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ.h>
namespace lhescheme {

class EvaluationKey {};

}

class LHE {
// A Base class for everyone 
private:
    long lambda;
    long tau;
public:
    LHE(long lambda, long tau);
};

