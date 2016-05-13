#ifndef POPULATION_INCLUDE
#define POPULATION_INCLUDE

#include <assert.h>
#include <cstddef>
#include <vector>
#include <limits>
#include <unordered_set>
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
        evaluated = new bool[npop_];

        for (int i=0;i<npop_;i++) {
            genes[i] = new T[ngene_];
            objectives[i] = new double[nobj_];
            constrains[i] = new double[nconst_];
            evaluated[i] = false;
        }
    }    
 
    ~Population() {
        for (int i=0;i<npop;i++) {
            delete [] genes[i];
            delete [] objectives[i];
            delete [] constrains[i];
        }
        delete [] genes;
        delete [] objectives;
        delete [] constrains;
        delete [] evaluated;
        
    }   
    
    void initialize();
    
    vector<vector<int>> nsga2(int generations) {
                
        //first evaluation
        printf("evaluating initial individuals...\n");
        //#pragma omp parallel for
        for (int i=0;i<npop;i++) {
            this->evaluator(genes[i],objectives[i],constrains[i]);
            printf("Initial ind[%d]: objectives=%f,%f: constrains=%f\n",i+1,objectives[i][0],objectives[i][1],constrains[i][0]);
        }
        printf("evaluating initial individuals...DONE\n");

        //for (int i=0;i<npop;i++) {
        //}

        exit(-19);

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
    
    void evolve(int generations, Population<T> &best) {
        //one generation more
        //1.- create offspring
        //printf("creating offspring\n");
        //printf("(1) Best individual at the moment: %f\n",bestIndividual->fitness);
        //Population offspring = Population(this);
        //printf("(1.0) Best individual at the moment: %f\n",bestIndividual->fitness);
        
        
        //first evaluation
        printf("evaluating initial individuals...%d\n",npop);
        //#pragma omp parallel for
        for (int i=0;i<npop;i++) {
            evaluate(i);
            //printf("evaluating[%d]: objectives=%f, constrains=%f\n",i,objectives[i][0],constrains[i][0]);
            //printAt(i);
        }
        printf("evaluating initial individuals...DONE\n");
        //Best
        copyTo(0,&best,0);
        for (int i=1;i<npop;i++) {
            if (bestTo(i,&best,0)) {
                copyTo(i,&best,0);
            }
        }
        printf("FIRST BEST: ");
        best.printAt(0);
        
       
        Population<T> offspring = Population<T>(npop*2,ngene,nobj,nconst,evaluator);

        double rnd;
        
        printf("Starting evolution\n");
        for (int g=0;g<generations;g++) {    
            //printf("Processing generation [%d]\n",g+1);
            //printf("Mating at generation [%d]\n",g+1);

            //copy current population half-size
            #pragma omp parallel for
            for (int i=0;i<npop;i++) {
                copyTo(i, &offspring, i);
            }
            
            //mate other half
            #pragma omp parallel for
            for (int i=npop;i<npop*2;i+=2) {
                //select 2 parents by Tournoument selection
                int p1 = select();
                int p2 = select();
                
                T* parent1 = genes[p1];
                T* parent2 = genes[p2];

                //crossover parents
                T* child1 = offspring.genes[i];
                T* child2 = offspring.genes[i+1];

                crossover(parent1,parent2,child1, child2);
                
                    //parent1->print();
                    //parent2->print();
                    //child1->print();
                    //child2->print();
                //offspring.evaluator(offspring.genes[npop+i],offspring.objectives[npop+i],offspring.constrains[npop+i]);
                //offspring.evaluator(offspring.genes[npop+i+1],offspring.objectives[npop+i+1],offspring.constrains[npop+i+1]);
                //this->printAt(p1,true);
                //this->printAt(p2,true);
                //offspring.printAt(npop+i,true);
                //offspring.printAt(npop+i+1,true);
                
                //both need to be evaluated
                offspring.evaluated[i] = false;
                offspring.evaluated[i+1] = false;
            }
            //printf("Mutating at generation [%d]\n",g+1);
            //2.- mutate
            #pragma omp parallel for
            for (int i=0;i<offspring.size();i++) {
                if (randomf() < mutationProb_) {
                    //offspring.get(i).print();
                    //printf("mutation prob=%f\n",mutationProb_);
                    //offspring.printAt(i);
                    mutate(offspring.genes[i]);
                    //offspring.printAt(i);

                    //mutated needs to be evaluated
                    offspring.evaluated[i] = false;

                    //offspring.get(i).print();
                    //printf("mutation prob=%f DONE!\n",indProb_);
                }
            }
            //3.- Evaluate
            //printf("Evaluating at generation [%d]\n",g+1);
            //printf("evaluation\n");
            //printf("(3) Best individual at the moment: %f\n",bestIndividual->fitness);
            //create array with 
            #pragma omp parallel for
            for (int i=0;i<offspring.size();i++) {
                //offspring.get(i).fitness = this->evaluator_(offspring.get(i));
                if (!offspring.evaluated[i]) {
                    offspring.evaluate(i);
                    //offspring.printAt(i);
                }
                //this->evaluator(offspring.genes[i],offspring.objectives[i],offspring.constrains[i]);
                //printf("evaluating[%d]: objectives=%f, constrains=%f\n",i,offspring.objectives[i][0],offspring.constrains[i][0]);
                //offspring.getRef(i).print();
            }
            //check best from offspring
            for (int i=0;i<offspring.size();i++) {
                if (offspring.bestTo(i,&best,0)) {
                    offspring.copyTo(i,&best,0);
                    printf("NEW BEST[gen:%d]: ",g+1);
                    best.printAt(0,false);
                }
            }
            //copy always the best --> elitism
            best.copyTo(0,this,0);

                //replace new generation by selection
            for (int i=1;i<npop;i++) {
                //select 2 parents by Tournoument selection
                int selected = offspring.select();
                
                offspring.copyTo(selected,this,i);
            }
            
            //printf("replacing new generation\n");
            double mean_fit = this->objectives[0][0];
            double min_fit = this->objectives[0][0];
            double max_fit = min_fit;

            double mean_cons = this->constrains[0][0];
            double min_cons = this->constrains[0][0];
            double max_cons = min_cons;
            
            for (int i=1;i<this->size();i++) {
                mean_fit += objectives[i][0];
                if (min_fit > objectives[i][0]) {
                    min_fit = objectives[i][0];
                }
                if (max_fit < objectives[i][0]) {
                    max_fit = objectives[i][0];
                }

                mean_cons += constrains[i][0];
                if (min_cons > constrains[i][0]) {
                    max_cons = constrains[i][0];
                }
                if (max_cons < constrains[i][0]) {
                    max_cons = constrains[i][0];
                }

            }
            mean_fit /= npop;
            mean_cons /= npop;
            
            printf("FITNESS[%d] Min=%f | Mean=%f | Max=%f. CONSTRAIN Min=%f | Mean=%f | Max=%f\n",g,min_fit,mean_fit,max_fit,min_cons,mean_cons,max_cons);
        }
    }
    


    void setup(int minvalue, int maxvalue,float crossoverProb = 0.8,float mutationProb = 0.2,float indpb = 0.01,float meta = 15,float tournamentSize = 5) {
        this->minvalue_ = minvalue;
        this->maxvalue_ = maxvalue;
        this->crossoverProb_ = crossoverProb;
        this->mutationProb_ = mutationProb;
        this->indProb_ = indpb;
        this->meta_ = meta;
        this->tournamentSize_ = tournamentSize;
    }
    
    size_t size() { return npop; }
    
    T **genes;
    double **objectives;
    double **constrains;

    size_t *rank;
    double *crowdDistance;
    bool *evaluated;

    
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

    void printAt(int index, bool printGene = false) {
        printIndividual(genes[index],ngene,objectives[index],nobj,constrains[index],nconst,0,0,index,printGene);
    }

    
    static void printIndividual(const T *gene,int ngene,const double *objectives,int nobj,const double *constrains,int nconst,int rank,double crowdDistance, int tag = 0, bool printGene = false) {
        if (tag>=0) {
            printf("%d,",tag);
        }
        for (int i=0;i<nobj;i++) {
            printf("%f,",i,objectives[i]);
        }
        for (int i=0;i<nconst;i++) {
            printf("%f,",i,constrains[i]);
        }
        printf(",%d,%f",rank,crowdDistance);

        if (printGene) {
            int acc = 0;
            printf("|GENE|");
            for (int i=0;i<ngene;i++) {
                printf("%d|",gene[i]);
                acc += gene[i];
            }
            printf("|mean gene=[%d]",acc/ngene);
        }
        
        printf("\n");
    }
    
    void randomize() {
        #pragma omp parallel for
        for (int i=0;i<npop;i++) {
            randomizeGene(genes[i],minvalue_,maxvalue_);
        }    
    }

    void initFromPopulation(imat &population,float initialPopulationProportion) {
        //imat populations has in each row an individual, plain row_order
        int tocopy = int(population.n_rows *initialPopulationProportion);
        tocopy = (tocopy > npop)?npop:tocopy;
        printf("population size:%d, external population size:%d, tocopy:%d, gensize=%d\n",npop,population.n_rows,tocopy,ngene);
        
        for (int i=0;i<tocopy;i++) {
            for (int j=0;j<ngene;j++) {
                genes[i][j] = population(i,j);
            }    
        }    

        for (int i=tocopy;i<npop;i++) {
            randomizeGene(genes[i],minvalue_,maxvalue_);
        }    

        printf("minvalue %d, maxvalue %d\n",minvalue_,maxvalue_);
        for (int i=0;i<npop;i++) {
            for (int j=0;j<ngene;j++) {
                assert(genes[i][j] >= minvalue_);
                assert(genes[i][j] <= maxvalue_);
            }    
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
    float meta_;
    int tournamentSize_;
    int minvalue_;
    int maxvalue_;

    int selectNSGA2();
    bool bestTo(int sel, Population *pop2, int selected);


    double (*evaluator)(T *,double *objs, double *consts);
    
    void randomizeGene(T *gene,int minvalue,int maxvalue) {
        //printf("randomize - IntArrayIndividual\n");
        for (int i=0;i<ngene;i++) {
            gene[i] = randomint(minvalue,maxvalue);
            //printf("gene[%d]=%d\n",i,gene[i]);
        }       
    }    
    
    void copyTo(int i, Population *other, int j);
    
    void evaluate(int i) {
        evaluator(this->genes[i],this->objectives[i],this->constrains[i]);
        this->evaluated[i] = true;
    }

};


//using IntPopulation = Population<int>;

template class Population<int>;


#endif
