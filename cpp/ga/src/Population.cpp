#include <stdio.h>
#include "Population.h"
#include "Random.h"
#include "Evaluator.h"
#include "Crossover.h"
#include "Mutation.h"
#include "Selection.h"

using namespace std;

Population::Population(size_t size, Individual& sample, evaluation_function_type evaluator): 
        sample_(sample),
        evaluator_(evaluator),
        bestIndividual(NULL) {

    individuals.clear();
    size_ = size;
    for (int i=0;i<size;i++) {
        Individual *newind = sample.clone();
        individuals.push_back(newind);
    }
}

Population::Population(Population *p): 
        sample_(p->sample_),
        evaluator_(p->evaluator_),
        bestIndividual(NULL) {

    individuals.clear();
    size_ = p->size();
    for (int i=0;i<size_;i++) {
        Individual *newind = p->sample_.clone();
        newind->copy(p->getPtr(i));
        individuals.push_back(newind);
    }
}

        
void Population::initialize() {
    randomize();
            
    //first evaluation
    #pragma omp parallel for
    for (int i=0;i<size_;i++) {
        individuals[i]->fitness = this->evaluator_(individuals[i]);
    }

    //printf("Initial population\n");
    //set best individual
    setNewBest(individuals[0]);
    for (int i=1;i<size_;i++) {
        //individuals[i]->print();
        if (bestIndividual->fitness < individuals[i]->fitness) {
            setNewBest(individuals[i]);
        }
    }

}

void Population::randomize() {
    vector<Individual *>::const_iterator it;
    #pragma omp parallel for private(it)
    for (it = individuals.begin() ; it < individuals.end(); ++it) {
        (*it)->randomize(minvalue_,maxvalue_);
    }    
}


Population::~Population() {
    vector<Individual *>::const_iterator it;
    #pragma omp parallel for private(it)
    for (it = individuals.begin() ; it < individuals.end(); ++it) {
        delete (*it);
    }
    
    if (bestIndividual != NULL) {
        printf("deleting best individual\n");
        //delete bestIndividual; 
    }
}


void Population::evolve() {
    //one generation more
    //1.- create offspring
    //printf("creating offspring\n");
    //printf("(1) Best individual at the moment: %f\n",bestIndividual->fitness);
    Population offspring = Population(this);
    //printf("(1.0) Best individual at the moment: %f\n",bestIndividual->fitness);
    
    double rnd;
    for (int i=0;i<size();i+=2) {
        //select 2 parents by Tournoument selection
        Individual* parent1 = select();
        Individual* parent2 = select();
        //crossover parents
        Individual* child1 = offspring.getPtr(i);
        Individual* child2 = offspring.getPtr(i+1);
        
        rnd = randomf();
        if (rnd < crossoverProb_) {
            //printf("performing crossover: [%f] < [%f]\n",rnd,crossoverProb_);
            crossover(parent1,parent2,child1, child2);
            //parent1->print();
            //parent2->print();
            //child1->print();
            //child2->print();

        } else {
            //printf("no crossover: [%f] < [%f]\n",rnd,crossoverProb_);
            child1->copy(parent1);
            child2->copy(parent2);
        }
        //printf("(1.1) Best individual at the moment: %f\n",bestIndividual->fitness);
    }
    //2.- mutate
    #pragma omp parallel for
    for (int i=0;i<size_;i++) {
        if (randomf() < mutationProb_) {
            //printf("mutation prob=%f\n",indProb_);
            //offspring.get(i).print();
            mutate(offspring.getPtr(i));
            //offspring.get(i).print();
        }
    }
    //3.- Evaluate
    //printf("evaluation\n");
    //printf("(3) Best individual at the moment: %f\n",bestIndividual->fitness);
    #pragma omp parallel for
    for (int i=0;i<size_;i++) {
        //offspring.get(i).fitness = this->evaluator_(offspring.get(i));
        this->evaluator_(offspring.getPtr(i));
        //offspring.getRef(i).print();
    }

    
    //replace new generation
    //printf("replacing new generation\n");
    double mean_fit = 0;
    double min_fit = offspring.getPtr(0)->fitness;
    double max_fit = min_fit;
    
    for (int i=0;i<size_;i++) {
        individuals[i]->copy(offspring.getPtr(i));
        if (bestIndividual->fitness < individuals[i]->fitness) {
            setNewBest(individuals[i]);
        }
        mean_fit += individuals[i]->fitness;
        if (min_fit > individuals[i]->fitness) {
            min_fit = individuals[i]->fitness;
        }
        if (max_fit < individuals[i]->fitness) {
            max_fit = individuals[i]->fitness;
        }
    }
    mean_fit /= size_;
    
    printf("FITNESS Min=%f | Mean=%f | Max=%f\n",min_fit,mean_fit,max_fit);
    
}

void Population::setNewBest(Individual* newbest) {
    if (bestIndividual == NULL) {
        bestIndividual = newbest->clone();
    }
    double prev = bestIndividual->fitness;
    bestIndividual->copy(newbest);
    
    printf("best new fitness=%f (was %f)\n",bestIndividual->fitness,prev);
    //bestIndividual->print();
    
}

Individual *Population::select() {
    return TournamentSelection(this,5);
    
    //int index1 = randomint(0,size_-1);
    //int index2 = randomint(0,size_-1);
    //if (individuals[index1]->fitness > individuals[index2]->fitness) {
        //return get(index1);
    //} else {
        //return get(index2);
    //}
}

Individual& Population::getRef(int i) {
    return *individuals[i];
}

Individual* Population::getPtr(int i) {
    return individuals[i];
}

void Population::crossover(Individual* parent1,Individual* parent2,Individual* child1, Individual* child2) {
    CrossoverUniform(parent1,parent2,child1,child2,0.5);
}

void Population::mutate(Individual* individual) {
    //GaussianMutation(individual,indProb_,minvalue_,maxvalue_,0,maxvalue_-minvalue_);
    SBXMutation<int>(individual,indProb_,10,minvalue_,maxvalue_);
}


void Population::setup(int minvalue, int maxvalue,float crossoverProb,float mutationProb,float indProb) {
    this->minvalue_ = minvalue;
    this->maxvalue_ = maxvalue;
    this->crossoverProb_ = crossoverProb;
    this->mutationProb_ = mutationProb;
    this->indProb_ = indProb;
}
