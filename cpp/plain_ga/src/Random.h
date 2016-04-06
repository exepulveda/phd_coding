#ifndef RANDOM_INCLUDE
#define RANDOM_INCLUDE

#include <cstdlib>
#include <random>

extern std::default_random_engine generator;
extern std::normal_distribution<double> norm;
extern std::uniform_real_distribution<double> uniform;

float randomf();
float randomf_gauss(float mu = 0.0, float sigma=1.0);
int randomint(int lower, int upper);

#endif
