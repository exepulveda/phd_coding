#include <cstdlib>

#include "Individual.h"
#include "Population.h"
#include "Random.h"

/* Routine for tournament selection, it creates a new_pop from old_pop by performing tournament selection and the crossover */
Individual* TournamentSelection(Population *population,int k) {
    int popSize = population->size();
    
    //select first
    int sel = randomint(0,popSize-1);
    Individual* selected = population->getPtr(sel);
    double max_fitness = selected->fitness;

    //printf("Selected: %d:%f\n",0,selected.fitness);

    
    for (int i=1;i<k;i++) {
        sel = randomint(0,popSize-1);
        Individual* ind = population->getPtr(sel);
        
        //printf("Select: %d:%f\n",i,ind.fitness);
        if (ind->fitness > max_fitness) {
            max_fitness = ind->fitness;
            selected = ind;
        }
    }

    //printf("Selected: %f\n",selected.fitness);
    
    return selected;
    
}
