#include <cstdlib>
#include <random>
#include "Random.h"

std::default_random_engine generator;
std::normal_distribution<double> norm(0.0,1.0);
std::uniform_real_distribution<double> uniform(0,1);


float randomf() { 
    return uniform(generator);
}

float randomf_gauss(float mu, float sigma) { 
    return norm(generator)*sigma + mu;
}

int randomint(int lower, int upper) { 
    return round(randomf() * (upper - lower)) + lower;
}

