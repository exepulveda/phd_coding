/* NSGA-II routine (implementation of the 'main' function) */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <unistd.h>
#include <armadillo>
#include <time.h>

#include "BlockCavingProblem.h"
#include "IntArrayIndividual.h"
#include "Population.h"

using namespace arma;

BlockCavingProblem bcp;

void individual2mat(IntArrayIndividual *individual, imat &schedule) {
    int k=0;
    for (int i=0;i<bcp.ndp;i++) {
        for (int j=0;j<bcp.nperiods;j++) {
            schedule(i,j) = individual->getInt(k);
            k++;
        }
    }
}


/* Routine to evaluate objective function values and constraints for an individual */
double evaluate(Individual *ind) {
    IntArrayIndividual *individual = (IntArrayIndividual *)ind;
    double obj[2];
    double cons[1];
    
    imat schedule(bcp.ndp,bcp.nperiods);

    individual2mat(individual,schedule);
    bcp.npv_cvar_npv(schedule,obj,cons);

    ind->fitness[0] = obj[0];
    ind->fitness[1] = obj[1];
    ind->constrains_[0] = cons[0];   
}


int main (int argc, char **argv) {
    int generations;
    int popsize;
    int seed;

    seed = 1634120;
    
    //loading blockcaving problem
    string config;

    if (argc >= 2) {
        popsize = atoi(argv[1]);
    } else {
        popsize = 100;
    }
    if (argc >= 3) {
        generations = atoi(argv[2]);
    } else {
        generations = 10;
    }
    if (argc >= 4) {
        config = string(argv[3]);
    } else {
        config = "/home/esepulveda/phd_coding/python/optimisation/pc1s1_B.json";
    }
    
    bcp.load(config);
    
    //for (int i=0;i<20;i++) {
        //printf("[%d] density=%f, nsr_average=%f\n",i,bcp.tonnage[i], bcp.nsr_average[i]);
    //}
    //for (int i=0;i<10;i++) {
        //for (int j=0;j<10;j++) {
            //printf("DP[%d] block[%d]=%d\n",i,j,bcp.getDrawPoint(i)->stream[j]);
        //}
    //}
    
    int numt = omp_get_max_threads();
    
    printf("Using %d threads\n",numt);
    
    const int ndp = 231;
    const int nperiods = 12;

    IntArrayIndividual sample(ndp*nperiods,2,1);
       
    Population population(popsize,sample,&evaluate);
    
    population.setup(0,200,0.9,0.2,0.05);

    population.nsga2(generations);

    printf("Best individual:\n");
    population.getBest().print();
    
    imat schedule(ndp,nperiods);
    individual2mat((IntArrayIndividual *)population.getBestPtr(),schedule);
    
    schedule.save("solution.mat",csv_ascii);
    
}
