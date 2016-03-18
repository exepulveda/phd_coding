#include "Individual.h"
#include "IntArrayIndividual.h"

IntArrayIndividual::IntArrayIndividual(const size_t N,int nobj,int nconst): Individual(N,nobj,nconst) { 
    this->gene = new int[N]; 
}

IntArrayIndividual::~IntArrayIndividual() {
    delete [] this->gene; 
}

void IntArrayIndividual::copy(Individual &other) {
    copy(&other);    
}

void IntArrayIndividual::copy(Individual *other) {
    Individual::copy(other);

    IntArrayIndividual *_o = (IntArrayIndividual *)other;
    for (int i=0;i<size();i++) {
        this->gene[i] = _o->gene[i];
    }    
}
