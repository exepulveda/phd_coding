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

imat individual2mat(IntArrayIndividual *individual) {
    imat schedule(bcp.ndp,bcp.nperiods);
    int k=0;
    for (int i=0;i<bcp.ndp;i++) {
        for (int j=0;j<bcp.nperiods;j++) {
            schedule(i,j) = individual->getInt(k);
            k++;
        }
    }
    return schedule;
}


/* Routine to evaluate objective function values and constraints for an individual */
double evaluate(Individual *ind) {
    IntArrayIndividual *individual = (IntArrayIndividual *)ind;
    double nsr;
    double deviation;
    double constrain;
    
    if (bcp.nsim > 0) {
        bcp.average_npv_tonnage_deviation(individual->gene,nsr,deviation,constrain);
    } else {
        imat schedule(bcp.ndp,bcp.nperiods);
        int k=0;
        for (int i=0;i<bcp.ndp;i++) {
            for (int j=0;j<bcp.nperiods;j++) {
                schedule(i,j) = individual->getInt(k);
                k++;
            }
        }
        bcp.single_npv_tonnage_deviation(schedule,nsr,constrain);
    }

    
    if (deviation < 0) {
        ind->fitness[0] = constrain;
    } else {
        ind->fitness[0] = nsr;
    }
    
    //create umatrix with blocks
    //printf("evaluating npv_cvar_npv [%d,%d]\n",bcp.ndp,bcp.nperiods);
    //time_t start,end;
    //time (&start);
    //bcp.npv_cvar_npv(ind->gene,ind->obj, ind->constr);

    //time (&end);

    //double elapsed_secs = difftime (end,start);
    
    //printf("evaluation obj1=%f, obj2=%f,const1=%f\n",ind->obj[0],ind->obj[1],ind->constr[0]);


    //ind->obj[0] = -ind->obj[0];
    //ind->obj[1] = -ind->obj[1];

    //ind->constr_violation = ind->constr[0];

    return ind->fitness[0];

    
    //imat schedule(bcp.ndp,bcp.nperiods);
    //int k=0;
    //for (int i=0;i<bcp.ndp;i++) {
    //    for (int j=0;j<bcp.nperiods;j++) {
    //        schedule(i,j) = ind->xreal[k++];
    //    }
    //}
    //bcp.npv_cvar_npv(schedule,ind->obj, ind->constr);

    


    //printf("obj1=[%d], obj2=[%f]\n",-ind->obj[0],ind->obj[1]);

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
        config = "/data/Linux/Documents/projects/newcrest/optimisation/data/pc1s1.json";
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

    IntArrayIndividual sample(ndp*nperiods);
    
    
    Population population(popsize,sample,&evaluate);
    
    population.setup(0,200,0.9,0.2,0.05);
    population.initialize();
    printf("Initial Best individual fitness=%f\n",population.getBest().fitness[0]);

    for (int i=0;i<generations;i++) {
        printf("Processing generation %d\n",i+1);
        population.evolve();
        printf("Best individual fitness=%f\n",population.getBest().fitness[0]);
    }

    printf("Best individual:\n");
    population.getBest().print();
    
    imat schedule = individual2mat((IntArrayIndividual *)population.getBestPtr());
    
    schedule.save("solution.mat",csv_ascii);
    
}
