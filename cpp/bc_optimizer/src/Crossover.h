#ifndef CROSSOVER_INCLUDE
#define CROSSOVER_INCLUDE

#include <cstdlib>

#include "Individual.h"

void CrossoverUniform(Individual* parent1,Individual* parent2,Individual* child1, Individual* child2,float indpb=0.5);


#endif
