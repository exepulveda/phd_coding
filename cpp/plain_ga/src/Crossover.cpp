#include "Crossover.h"
#include "Individual.h"
#include "Random.h"

void CrossoverUniform(Individual* parent1,Individual* parent2,Individual* child1, Individual* child2,float indpb) {
    int size1 = parent1->size();
    int size2 = parent2->size();
    
    int size = (size1>size2)?size2:size1;
    
    for (int i=0;i<size;i++) {
        if (randomf() < indpb) {
            child1->set(i, parent2->getInt(i));
            child2->set(i, parent1->getInt(i));
        } else {
            child1->set(i, parent1->getInt(i));
            child2->set(i, parent2->getInt(i));
        }
    }
}
