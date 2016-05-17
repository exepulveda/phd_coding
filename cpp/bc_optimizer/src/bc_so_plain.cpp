#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <armadillo>
#include <time.h>
#include <sys/stat.h>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "BlockCavingProblem.h"
#include "Population.hpp"

using namespace arma;
using namespace std;

static BlockCavingProblem bcp;
static int currentSimulation = -1;

template<typename T>
void individual2mat(T *individual, imat &schedule) {
    //row order
    int k=0;
    for (int i=0;i<bcp.ndp;i++) {
        for (int j=0;j<bcp.nperiods;j++) {
            schedule(i,j) = individual[k];
            k++;
        }
    }
}

template<typename T>
void evaluateNSR(T *individual,double *objs, double *consts) {
    double nsr;
    double constrain;
    
    if (currentSimulation >= 0) {
        bcp.single_npv_tonnage_deviation(individual,nsr,constrain,currentSimulation);
    } else {
        bcp.average_npv_tonnage(individual,nsr,constrain);
    }

    //printf("nsr=%f,constrain=%f\n",nsr,constrain);

    objs[0] = nsr; 
    consts[0] = constrain;

}

namespace pt = boost::property_tree;

using namespace boost::property_tree;
using namespace pt::json_parser;
using namespace boost::filesystem;

inline bool existFile(const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

int main (int argc, char **argv) {
    int caseEval = 1;
    char filename[1000];
    //loading blockcaving problem
    string config;
    string gaConfig;
    string destination;
    // Declare the supported options.

    printf("argc=%d\n",argc);
    for (int i=0;i<argc;i++) {
        printf("argv[%d]=[%s]\n",i,argv[i]);
    }


    if (argc < 3) {
        printf("use: %s config destination [gaconfig] [sim] [case]\n",argv[0]);
        return (-1);
    }
    
    //check if files exist
    config = string(argv[1]);
    
    if (existFile(config)) {
        printf("Using config file [%s]\n",config.c_str());
    } else {
        printf("Config file [%s] does not exist\n",config.c_str());
        return(-4);
    }

    int popsize = 100;
    int generations = 50;
    int seed = 1634120;
    float crossover = 1.0;
    float mutation = 0.2;
    float mutationGenePbl = 0.05;
    float mutationEta = 15;
    int tournamentSize = 5;

    if (argc > 3) {
        gaConfig = string(argv[3]);
        printf("Using gaConfig=%s\n",gaConfig.c_str());
        if (existFile(gaConfig)) {
            ptree pt;
            json_parser::read_json(gaConfig, pt);
            
            popsize = pt.get<int>("population");
            generations = pt.get<int>("generations");
            seed = pt.get<int>("seed");
            crossover = pt.get<float>("crossover.probability");
            mutation = pt.get<float>("mutation.probability");
            mutationGenePbl = pt.get<float>("mutation.geneProbability");
            mutationEta = pt.get<float>("mutation.sbx_type.eta");
            tournamentSize = pt.get<int>("tournament_size");
        } else {
            printf("GA config file does not exist\n");
            return(-3);
        }
    }
        
    if (argc > 4) {
        currentSimulation = atoi(argv[4]);
    } else {
        currentSimulation = -1;
    }

    if (argc > 5) {
        caseEval = atoi(argv[5]);
    }


    //delete destination folder
    destination = string(argv[2]);
    printf("destination=%s...\n",argv[2]);

    printf("loading BCP problem config...\n");
    bcp.load(config);
    printf("loading BCP problem config...DONE!\n");
    
    printf("Max threads %d\n",omp_get_max_threads());
    printf("Using %d threads\n",omp_get_num_procs());
    
    const int ndp = bcp.ndp;
    const int nperiods = bcp.nperiods;
    
    void (*evaluateFunction)(int *gene,double *objs, double *consts);
    
    switch (caseEval) {
        case 1:
            evaluateFunction = &evaluateNSR;
            break;
        case 2:
            evaluateFunction = &evaluateNSR;
            break;
        case 3:
            evaluateFunction = &evaluateNSR;
            break;
        default:
            evaluateFunction = &evaluateNSR;
            break;
    }

    Population<int> population(popsize,ndp*nperiods,1,1,evaluateFunction);
    
    population.setup(0,bcp.maxExtraction,1.0,mutation,mutationGenePbl,mutationEta,tournamentSize);

    printf("Initial population from random\n");
    population.randomize();

    //best individual
    Population<int> bestIndividual(1,ndp*nperiods,1,1,evaluateFunction);

    //evolution
    fflush(stdout);
    
    population.evolve(generations,bestIndividual);
    
    fflush(stdout);    
    
    printf("Best individual:\n");
    imat schedule(ndp,nperiods);    
    
    individual2mat(bestIndividual.genes[0],schedule);
    
    
    sprintf(filename,"%s-solution.mat",destination.c_str());
    schedule.save(filename,csv_ascii);
    
    printf("Saving population...\n");
    sprintf(filename,"%s-popultation.csv",destination.c_str());
    FILE *fd = fopen(filename,"w");
    //the best always is copied
    for (int j=0;j<ndp*nperiods;j++) {
        if (j ==0) {
            fprintf(fd,"%d",bestIndividual.genes[0][j]);
        } else {
            fprintf(fd,",%d",bestIndividual.genes[0][j]);
        }
    }
    fprintf(fd,"\n");
    
    for (int i=0;i<popsize;i++) {
        for (int j=0;j<ndp*nperiods;j++) {
            if (j ==0) {
                fprintf(fd,"%d",population.genes[i][j]);
            } else {
                fprintf(fd,",%d",population.genes[i][j]);
            }
        }
        fprintf(fd,"\n");
    }
    
    fclose(fd);
}
