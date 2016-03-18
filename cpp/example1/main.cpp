#include <omp.h>

#include "Individual.h"
#include "IntArrayIndividual.h"
#include "Population.h"

double evaluateOnes(Individual *ind) {
    IntArrayIndividual *individual = (IntArrayIndividual *)ind;
    
    int ones = 0;
    for (int i=0;i<individual->size();i++) {
        //printf("Ind[%d]=%d\n",i,individual->getInt(i));
        ones += individual->getInt(i);
    }
    
    individual->fitness[0] = ones;
    return ones;
}


int main() {
    int numt = omp_get_max_threads();
    
    printf("Using %d threads\n",numt);
    
    const int ndp = 231;
    const int nperiods = 12;
    IntArrayIndividual sample(20);
    
    int generations = 1000;
    
    Population population(50,sample,&evaluateOnes);
    
    population.setup(0,9,0.9,0.2,0.05);
    population.initialize();
    printf("Initial Best individual fitness=%f\n",population.getBest().fitness[0]);

    for (int i=0;i<generations;i++) {
        printf("Processing generation %d\n",i+1);
        population.evolve();
        //printf("Best individual fitness=%f\n",population.getBest().fitness);
    }

    printf("Best individual:\n");
    population.getBest().print();

}
