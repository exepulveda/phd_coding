/* NSGA-II routine (implementation of the 'main' function) */

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

BlockCavingProblem bcp;

template<typename T>
void evaluateNSRDeviation(T *individual,double *objs, double *consts) {
    double nsr;
    double dev;
    double constrain;
    
    bcp.average_nsr_tonnage_deviation(individual,nsr,dev,constrain);

    objs[0] = nsr;
    objs[1] = -dev; //minimize
    consts[0] = constrain;   
}

template<typename T>
void evaluateNSRVariance(T *individual,double *objs, double *consts) {
    double nsr;
    double var;
    double constrain;
    
    bcp.average_nsr_variance(individual,nsr,var,constrain);

    objs[0] = nsr;
    objs[1] = -var; //minimize
    consts[0] = constrain;   
}

template<typename T>
void evaluateNSRCVaR(T *individual,double *objs, double *consts) {
    double nsr;
    double cvar;
    double constrain;
    
    bcp.average_nsr_cvar(individual,nsr,cvar,constrain);

    objs[0] = nsr;
    objs[1] = cvar;
    consts[0] = constrain;   
}

template<typename T>
void evaluateNSRVaR(T *individual,double *objs, double *consts) {
    double nsr;
    double var;
    double constrain;
    
    bcp.average_nsr_var(individual,nsr,var,constrain);

    objs[0] = nsr;
    objs[1] = var;
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
        printf("use: %s config destination [gaconfig] [case]\n",argv[0]);
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
        caseEval = atoi(argv[4]);
    }

    string initialIndividual;
    if (argc > 5) {
        initialIndividual = string(argv[5]);
    }
    int initialIndividualCopies = int(popsize*0.1); //10%
    if (argc > 6) {
        initialIndividualCopies = atoi(argv[6]);
    } 

    //delete destination folder
    destination = string(argv[2]);
    printf("destination=%s...\n",argv[2]);
    //system::error_code ec;
    //remove_all(destination,ec); 
    //create an empty destinatino folder
    //if (!create_directories(destination)) {
    //    printf("The destination folder could not be created\n");
    //    return(-2);
    //}
    //printf("destination DONE!\n");

    printf("loading BCP problem config...\n");
    bcp.load(config);
    printf("loading BCP problem config...DONE!\n");
    
    //for (int i=0;i<20;i++) {
        //printf("[%d] density=%f, nsr_average=%f\n",i,bcp.tonnage[i], bcp.nsr_average[i]);
    //}
    //for (int i=0;i<10;i++) {
        //for (int j=0;j<10;j++) {
            //printf("DP[%d] block[%d]=%d\n",i,j,bcp.getDrawPoint(i)->stream[j]);
        //}
    //}
    
    printf("Max threads %d\n",omp_get_max_threads());
    printf("Using %d threads\n",omp_get_num_procs());
    
    const int ndp = bcp.ndp;
    const int nperiods = bcp.nperiods;
    
    void (*evaluateFunction)(int *gene,double *objs, double *consts) = &evaluateNSRVariance<int>;
    
    switch (caseEval) {
        case 1:
            evaluateFunction = &evaluateNSRDeviation;
            break;
        case 2:
            evaluateFunction = &evaluateNSRVariance;
            break;
        case 3:
            evaluateFunction = &evaluateNSRVaR;
            break;
        case 4:
            evaluateFunction = &evaluateNSRCVaR;
            break;
    }

    Population<int> population(popsize,ndp*nperiods,2,1,evaluateFunction);
    
    population.setup(0,bcp.maxExtraction,1.0,mutation,mutationGenePbl,mutationEta,tournamentSize);

    if (initialIndividual.length() == 0) {
        printf("Initial population from random\n");
        population.randomize();
    } else {
        printf("Initial individual from external solution\n");
        imat ipop;
        ipop.load(initialIndividual,csv_ascii);
        population.initFromPopulation(ipop,initialIndividualCopies);
    }

    vector<vector<int>> front = population.nsga2(generations);

    printf("Best individuals:\n");

    population.printFront(front[0]);
    
    imat schedule(ndp,nperiods);
    for (int i : front[0]) {
        sprintf(filename,"solution-mo-%d.mat",i);
        
        boost::filesystem::path pout(destination);
        pout /= filename;
        
        bcp.individual2mat(population.genes[i],schedule);
        
        schedule.save(pout.string(),csv_ascii);
        
        /*
        FILE *fs = fopen(pout.string().c_str(),"w");
        
        int k=0;
        for (int dp=0;dp<ndp;dp++) {
            for (int p=0;p<nperiods;p++) {
                if (p == 0) {
                    fprintf(fs,"%d",population.genes[i][k]);
                } else {
                    fprintf(fs,",%d",population.genes[i][k]);
                }
                k++;
            }
            fprintf(fs,"\n");
        }
        
        fclose(fs);
        */
    }
}
