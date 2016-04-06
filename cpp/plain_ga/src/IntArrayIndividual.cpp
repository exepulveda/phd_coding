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

void IntArrayIndividual::print(int tag, bool printGene) {
    if (tag>=0) {
        printf("Individual[%d]|",tag);
    }
    if (printGene) {
        for (int i=0;i<size();i++) {
            printf("%d|",gene[i]);
        }
    }
    printf("fitness|");
    for (int i=0;i<nobj;i++) {
        printf("[%d]=%f|",i,fitness[i]);
    }
    printf("constrains|");
    for (int i=0;i<nconst;i++) {
        printf("[%d]=%f|",i,constrains_[i]);
    }
    printf(". Rank=%d, CrowdDistance=%f\n",this->rank,this->crowdDistance);
}
