#ifndef POPULATION_INCLUDE
#define POPULATION_INCLUDE

#include <assert.h>
#include <cstddef>
#include <vector>
#include <limits>
#include <armadillo>
using namespace arma;

#include "Random.h"

using namespace std;

const static double inf = numeric_limits<double>::infinity();

template <typename T>
class Population {
public:
    Population(size_t npop_, size_t ngene_, size_t nobj_, size_t nconst_, double (*evaluator_)(T *gene,double *objs, double *consts)): 
            npop(npop_), 
            ngene(ngene_), 
            nobj(nobj_), 
            nconst(nconst_), 
            evaluator(evaluator_) 
            {

        //init memory
        genes = new T*[npop_];
        objectives = new double*[npop_];
        constrains = new double*[npop_];
        rank = new size_t[npop_];
        crowdDistance = new double[npop_];

        for (int i=0;i<npop_;i++) {
            genes[i] = new T[ngene_];
            objectives[i] = new double[nobj_];
            constrains[i] = new double[nconst_];
        }
    }    
 
    ~Population() {
        for (int i=0;i<npop;i++) {
            delete [] genes[i];
            delete [] objectives[i];
            delete [] constrains[i];
        }
        delete genes;
        delete objectives;
        delete constrains;
        
    }   
    
    void evolve();
    void initialize();
    
    vector<vector<int>> nsga2(int generations) {
                
        //first evaluation
        printf("evaluating initial individuals...\n");
        #pragma omp parallel for
        for (int i=0;i<npop;i++) {
            this->evaluator(genes[i],objectives[i],constrains[i]);
        }
        printf("evaluating initial individuals...DONE\n");

        //calculate first rank
        vector< vector < int > > fronts;

        printf("first calculateFrontRank...\n");
        calculateFrontRank(fronts);
        printf("first calculateFrontRank...DONE\n");

        //printAllFronts(fronts);

        T *parent1;
        T *parent2;
        T *child1;
        T *child2;
        
        printf("creating ofgspring of size=%d,gensize=%d\n",npop*2,ngene);

        Population<T> offspring = Population<T>(npop*2,ngene,nobj,nconst,evaluator);
        offspring.setup(minvalue_,maxvalue_,crossoverProb_,mutationProb_,indProb_);
                   
        size_t geneSizeBytes = sizeof(T)*ngene;
        size_t objSizeBytes = sizeof(double)*nobj;
        size_t constSizeBytes = sizeof(double)*nconst;
                   
        for (int g=0;g<generations;g++) {
            printf("processing generation %d/%d\n",g+1,generations);
            //update offspring

            
            #pragma omp parallel for private(parent1,parent2,child1,child2)
            for (int i=0;i<npop;i+=2) {
                //select 2 parents by Tournoument selection
                parent1 = genes[selectNSGA2()];
                parent2 = genes[selectNSGA2()];
                //crossover parents
                child1 = offspring.genes[i];
                child2 = offspring.genes[i+1];
                
                //printf("performing crossover\n");
                offspring.crossover(parent1,parent2,child1, child2);
                //parent1->print();
                //parent2->print();
                //child1->print();
                //child2->print();

                //2.- mutate
                //printf("performing mutation\n");
                if (randomf() < mutationProb_) {
                    //printf("mutation prob=%f\n",indProb_);
                    //offspring.get(i).print();
                    offspring.mutate(child1);
                    //offspring.get(i).print();
                }
                if (randomf() < mutationProb_) {
                    //printf("mutation prob=%f\n",indProb_);
                    //offspring.get(i).print();
                    offspring.mutate(child2);
                    //offspring.get(i).print();
                }
                //printf("evaluating\n");
                //3.- Evaluate
                offspring.evaluator(child1,offspring.objectives[i],offspring.constrains[i]);
                offspring.evaluator(child2,offspring.objectives[i+1],offspring.constrains[i+1]);

            }

            //add original population to offspring
            printf("adding original population to offspring [%d:%d][%d:%d:%d]\n",npop,offspring.size(),geneSizeBytes,objSizeBytes,constSizeBytes);
            #pragma omp parallel for
            for (int i=npop;i<offspring.size();i++) {
                //memcpy(offspring.genes[i],genes[i-npop],geneSizeBytes);
                //memcpy(&offspring.objectives[i][0],&objectives[i-npop][0],objSizeBytes);
                //memcpy(&offspring.constrains[i][0],&constrains[i-npop][0],constSizeBytes);
                
                copyTo(i-npop, &offspring, i);

                //printf("copying from population[%d]:\n",i-npop);
                //Population::printIndividual(genes[i-npop],ngene,objectives[i-npop],nobj,constrains[i-npop],nconst,rank[i-npop],crowdDistance[i-npop], i-npop,true);                
                //printf("copying to offspring[%d]:\n",i);
                //Population::printIndividual(offspring.genes[i],ngene,offspring.objectives[i],nobj,offspring.constrains[i],nconst,offspring.rank[i],offspring.crowdDistance[i], i,true);                
            }
            
            //purge duplicated
            printf("purging duplicated\n");
            offspring.purgeDuplicated();

            //printf("offspring:\n");
            //for (int i=0;i<offspring.size();i++) {
            //    Population::printIndividual(offspring.genes[i],ngene,offspring.objectives[i],nobj,offspring.constrains[i],nconst,offspring.rank[i],offspring.crowdDistance[i], i,false);                
            //}


            printf("refreshing rank\n");
            //calculate rank and crwodistance
            fronts.clear();
            printf("calculateFrontRank...\n");
            offspring.calculateFrontRank(fronts);
            printf("calculateFrontRank...DONE\n");

            //offspring.printAllFronts(fronts);


            fronts = selectFront(fronts,npop);

            //apply new front
            int i = 0;
            for(int f=0;f<fronts.size();f++) {
                vector<int> &front = fronts[f];
                for(int j=0;j<front.size();j++) {
                    //memcpy(genes[i],offspring.genes[front[j]],geneSizeBytes);
                    //memcpy(objectives[i],offspring.objectives[front[j]],objSizeBytes);
                    //memcpy(constrains[i],offspring.constrains[front[j]],constSizeBytes);
                    //crowdDistance[i] = offspring.crowdDistance[front[j]];

                    offspring.copyTo(front[j], this, i);
                    rank[i] = f;
                    front[j] = i;
                    i++;
                }
            }
            
            //printAllFronts(fronts);
            printFront(fronts[0],0);
            //printFront(fronts[1],1);

            assert (i == npop);        
        }
        
        return fronts;
    }


    void setup(int minvalue, int maxvalue,float crossoverProb = 0.8,float mutationProb = 0.2,float indpb = 0.01) {
        this->minvalue_ = minvalue;
        this->maxvalue_ = maxvalue;
        this->crossoverProb_ = crossoverProb;
        this->mutationProb_ = mutationProb;
        this->indProb_ = indpb;
    }
    
    size_t size() { return npop; }
    
    T **genes;
    double **objectives;
    double **constrains;

    size_t *rank;
    double *crowdDistance;

    
    size_t npop;
    size_t ngene;
    size_t nobj;
    size_t nconst;

    void printAllFronts(vector<vector<int>> &front) {
        int i=0;
        for (vector<int> &f : front) {
            printf("Front[%d], size=%d. Elements:\n",i,f.size());
            printFront(f,i,false);
            i++;
        }
    }
    void printFront(vector<int> &f,int tag = 0, bool printGene = false) {
        for (int p : f) {
            printIndividual(genes[p],ngene,objectives[p],nobj,constrains[p],nconst,rank[p],crowdDistance[p],p, printGene);
        }
    }    
    
    static void printIndividual(const T *gene,int ngene,const double *objectives,int nobj,const double *constrains,int nconst,int rank,double crowdDistance, int tag = 0, bool printGene = false) {
        if (tag>=0) {
            printf("%d,",tag);
        }
        if (printGene) {
            for (int i=0;i<ngene;i++) {
                printf("%d|",gene[i]);
            }
        }
        for (int i=0;i<nobj;i++) {
            printf("%f,",i,objectives[i]);
        }
        for (int i=0;i<nconst;i++) {
            printf("%f,",i,constrains[i]);
        }
        printf(",%d,%f\n",rank,crowdDistance);
    }
    
    void randomize() {
        #pragma omp parallel for
        for (int i=0;i<npop;i++) {
            randomizeGene(genes[i],minvalue_,maxvalue_);
        }    
    }

    void initFromPopulation(imat &population) {
        int tocopy = (population.n_rows > npop)?npop:population.n_rows;
        printf("population size:%d, external population size:%d, tocopy:%d\n",npop,population.n_rows,tocopy);
        for (int i=0;i<tocopy;i++) {
            for (int j=0;j<ngene;j++) {
                genes[i][j] = population(i,j);
            }    
        }    
        
        for (int i=tocopy;i<npop;i++) {
            randomizeGene(genes[i],minvalue_,maxvalue_);
        }    
    }
    
private:
    
    int select();
    void crossover(const T *parent1,const T *parent2,T *child1,T *child2);
    void mutate(T *individual);
    void purgeDuplicated();

    vector<vector<int>> selectFront(vector<vector<int>> &front,int newSize);
    void calculateFrontRank(vector<vector<int>> &front);
    void calculateCrowdingDistance(vector<int> front);

    static void SBXMutation(T *individual, int ngene, float indpb,float eta_m,T minvalue,T maxvalue);
    int tournamentSelection(int k);
    int tournamentSelectionNSGA2(int k);
    static void crossoverUniform(const T *parent1,const T *parent2,T *child1,T *child2,int ngene,float indpb);
    static void gaussianMutation(T *individual,int ngene, float indpb,T minvalue,T maxvalue, float mu, float sigma);
    bool dominates(int p,int q);

    
    float mutationProb_;
    float crossoverProb_;
    float indProb_;
    int minvalue_;
    int maxvalue_;

    int selectNSGA2();

    double (*evaluator)(T *,double *objs, double *consts);
    
    void randomizeGene(T *gene,int minvalue,int maxvalue) {
        //printf("randomize - IntArrayIndividual\n");
        for (int i=0;i<ngene;i++) {
            gene[i] = randomint(minvalue,maxvalue);
            //printf("gene[%d]=%d\n",i,gene[i]);
        }       
    }    
    
    void copyTo(int i, Population *other, int j);
};


//using IntPopulation = Population<int>;

template class Population<int>;


#endif
