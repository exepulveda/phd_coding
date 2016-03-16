#include "Individual.h"
#include "IntArrayIndividual.h"

IntArrayIndividual::IntArrayIndividual(const size_t N): Individual(N) { 
    this->gene = new int[N]; 
}

IntArrayIndividual::~IntArrayIndividual() {
    delete [] this->gene; 
}

void IntArrayIndividual::copy(Individual &other) {
    IntArrayIndividual *_o = (IntArrayIndividual *)&other;
    for (int i=0;i<size();i++) {
        this->gene[i] = _o->gene[i];
    }
    
    this->fitness = other.fitness;
    
}

void IntArrayIndividual::copy(Individual *other) {
    copy(*other);
}
