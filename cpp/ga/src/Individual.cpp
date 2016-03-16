#include "Individual.h"

void Individual::copy(Individual &other) { 
    fitness = other.fitness; 
    size_ = other.size_; 
}

void Individual::copy(Individual *other) { 
    fitness = other->fitness; 
    size_ = other->size_; 
}
