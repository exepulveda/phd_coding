/* NSGA-II routine (implementation of the 'main' function) */

# include <vector>
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <unistd.h>
#include <armadillo>
#include <time.h>
# include <boost/program_options.hpp>

#include "BlockCavingProblem.h"
#include "IntArrayIndividual.h"
#include "Population.h"
#include "FrontRank.h"

using namespace arma;
using namespace std;

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
    //bcp.npv_cvar_npv(individual->gene,obj,cons);


    ind->fitness[0] = obj[0];
    ind->fitness[1] = obj[1];
    ind->constrains_[0] = cons[0];   
}

double evaluateNSRDeviation(Individual *ind) {
    IntArrayIndividual *individual = (IntArrayIndividual *)ind;
    double nsr;
    double dev;
    double constrain;
    
    bcp.average_npv_tonnage_deviation(individual->gene,nsr,dev,constrain);

    ind->fitness[0] = nsr;
    ind->fitness[1] = -dev; //minimize
    ind->constrains_[0] = constrain;   
}


namespace po = boost::program_options;

int main (int argc, char **argv) {
    int generations;
    int popsize;
    int seed;
    int caseEval;
    char filename[1000];
    //loading blockcaving problem
    string config;

    // Declare the supported options.
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("npop", po::value<int>(&popsize)->default_value(100), "size of population")
        ("ngen", po::value<int>(&generations)->default_value(100), "number of generations")
        ("seed", po::value<int>(&seed)->default_value(1634120), "seed of random numbers")
        ("case", po::value<int>(&caseEval)->default_value(1), "case to run")
        ("config", po::value<string>(&config)->default_value("/home/esepulveda/phd_coding/python/optimisation/pc1s1_B.json"), "configuration file")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    

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
    
    
    evaluation_function_type evaluator;
    switch (caseEval) {
        case 1:
            evaluator = &evaluate;
            break;
        case 2:
            evaluator = &evaluateNSRDeviation;
            break;
        default:
            evaluator = &evaluate;
            break;
    }
    
    Population population(popsize,sample,evaluator);
    
    population.setup(0,200,1.0,0.2,0.05);

    vector< vector < int > > front = population.nsga2(generations);

    printf("Best individuals:\n");

    printFront(&population,front[0]);
    
    for (int i : front[0]) {
        imat schedule(ndp,nperiods);
        individual2mat((IntArrayIndividual *)population.getPtr(i),schedule);
    
        sprintf(filename,"solution-mo-%d.mat",i);
        schedule.save(string(filename),csv_ascii);
    }
}
