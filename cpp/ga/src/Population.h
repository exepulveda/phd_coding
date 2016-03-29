#ifndef POPULATION_INCLUDE
#define POPULATION_INCLUDE

#include <cstddef>
#include <vector>

#include "Individual.h"

using namespace std;

typedef double (*evaluation_function_type)(Individual *);

class Population {
public:
    Population(size_t size, Individual& sample, evaluation_function_type evaluator);
    Population(Population *p);
    
    ~Population();
    
    void evolve();
    void initialize();
    
    vector< vector < int > > nsga2(int ngen);

    
    Individual& getRef(int i);
    Individual* getPtr(int i);
    
    void setup(int minvalue, int maxvalue,float crossoverProb = 0.8,float mutationProb = 0.2,float indpb = 0.01);
        
    size_t size() { return size_; }
    
    Individual& getBest() { return *bestIndividual; }
    Individual* getBestPtr() { return bestIndividual; }

    //void set(int i, Individual& ind) { individuals[i] = ind; }
    
private:
    void randomize();
    vector<Individual *> individuals;
    Individual& sample_;
    
    Individual* select();
    void crossover(Individual* parent1,Individual* parent2,Individual* child1, Individual* child2);
    void mutate(Individual* individual);
    void append(Population *other);
    Population *purgeDuplicated();
    float mutationProb_;
    float crossoverProb_;
    float indProb_;
    int minvalue_;
    int maxvalue_;

    Individual* selectNSGA2();

    
    evaluation_function_type evaluator_;
    
    Individual* bestIndividual;

    void setNewBest(Individual* newbest);
    
    size_t size_;

};

#endif
