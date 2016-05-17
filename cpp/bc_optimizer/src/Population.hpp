#ifndef POPULATION_INCLUDE
#define POPULATION_INCLUDE

#include <assert.h>
#include <cstddef>
#include <vector>
#include <limits>
#include <unordered_set>

#define ARMA_NO_DEBUG

#include <armadillo>

#include <stdio.h>
#include <algorithm>
#include <cstring>

#include "Random.h"

using namespace arma;

using namespace std;

const static double inf = numeric_limits<double>::infinity();

bool myfunction(pair<double,int> i,pair<double,int> j) { return (i.first<j.first); }


template <typename T>
class Population {
public:
    Population(size_t npop_, size_t ngene_, size_t nobj_, size_t nconst_, void (*evaluator_)(T *gene,double *objs, double *consts)): 
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
    
    vector<vector<int>> nsga2(int generations) {
                
        //first evaluation
        printf("evaluating initial individuals...\n");
        #pragma omp parallel for
        for (int i=0;i<npop;i++) {
            this->evaluator(genes[i],objectives[i],constrains[i]);
            printf("Initial ind[%d]: objectives=%f,%f: constrains=%f\n",i+1,objectives[i][0],objectives[i][1],constrains[i][0]);
        }
        printf("evaluating initial individuals...DONE\n");
        fflush(stdout);
        //for (int i=0;i<npop;i++) {
        //}

        //calculate first rank
        vector< vector < int > > fronts;

        printf("first calculateFrontRank...\n");
        calculateFrontRank(fronts);
        printf("first calculateFrontRank...DONE\n");
        fflush(stdout);

        //printAllFronts(fronts);

        T *parent1;
        T *parent2;
        T *child1;
        T *child2;
        
        printf("creating offspring of size=%d,gensize=%d\n",npop*2,ngene);

        Population<T> offspring = Population<T>(npop*2,ngene,nobj,nconst,evaluator);
        offspring.setup(minvalue_,maxvalue_,crossoverProb_,mutationProb_,indProb_,meta_,tournamentSize_);
                   
        for (int g=0;g<generations;g++) {
            printf("processing generation %d/%d\n",g+1,generations);
            fflush(stdout);
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
            //printf("adding original population to offspring [%d:%d][%d:%d:%d]\n",npop,offspring.size(),geneSizeBytes,objSizeBytes,constSizeBytes);
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


            //printf("refreshing rank\n");
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
            fflush(stdout);

            //assert (i == npop);        
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
        #pragma omp parallel for
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
        
        fflush(stdout);       
       
        Population<T> offspring = Population<T>(npop*2,ngene,nobj,nconst,evaluator);

        printf("Starting evolution\n");
        fflush(stdout);
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
                    min_cons = constrains[i][0];
                }
                if (max_cons < constrains[i][0]) {
                    max_cons = constrains[i][0];
                }

            }
            mean_fit /= npop;
            mean_cons /= npop;
            
            printf("FITNESS[%d] Min=%f | Mean=%f | Max=%f. CONSTRAIN Min=%f | Mean=%f | Max=%f\n",g,min_fit,mean_fit,max_fit,min_cons,mean_cons,max_cons);
            fflush(stdout);
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

    
    void printIndividual(const T *gene,int ngene,const double *objectives,int nobj,const double *constrains,int nconst,int rank,double crowdDistance, int tag = 0, bool printGene = false) {
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
    /*
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
    bool bestTo(int sel, Population *pop2, int selected);
    void copyTo(int i, Population *other, int j);
    */
    
    float mutationProb_;
    float crossoverProb_;
    float indProb_;
    float meta_;
    int tournamentSize_;
    int minvalue_;
    int maxvalue_;
    void (*evaluator)(T *,double *objs, double *consts);
    
    void randomizeGene(T *gene,int minvalue,int maxvalue) {
        //printf("randomize - IntArrayIndividual\n");
        for (int i=0;i<ngene;i++) {
            gene[i] = randomint(minvalue,maxvalue);
            //printf("gene[%d]=%d\n",i,gene[i]);
        }       
    }    
    
    
    void evaluate(int i) {
        evaluator(this->genes[i],this->objectives[i],this->constrains[i]);
        this->evaluated[i] = true;
    }

    int selectNSGA2() {
        return tournamentSelectionNSGA2(tournamentSize_);
    };

    void initialize() {
        randomize();
                
        //first evaluation
        #pragma omp parallel for
        for (int i=0;i<npop;i++) {
            this->evaluator(genes[i],objectives[i],constrains[i]);
        }
    }

    int select() {
        return tournamentSelection(tournamentSize_);
    }

    void crossover(const T *parent1,const T *parent2,T *child1,T *child2) {
        crossoverUniform(parent1,parent2,child1,child2,ngene,0.5);
    }

    void mutate(T *individual) {
        SBXMutation(individual,ngene,indProb_,meta_,minvalue_,maxvalue_);
    }

    void purgeDuplicated() {
        //purging is a trick to give an inf rank to duplicates
        //create a set
        T *pi;
        T *pj;
        int nduplicated = 0;

        double minConst = constrains[0][0];
        for (int i=1;i<size();i++) {
            if (constrains[i][0] < minConst) {
                minConst = constrains[i][0];
            }
        }
        
        
        for (int i=0;i<size();i++) {
            pi = genes[i];
            bool unique = true;
            for (int j=i+1;j<size();j++) {
                pj = genes[j];
                if (memcmp(pi,pj,sizeof(T)*ngene) == 0) {
                    unique = false;
                    break;
                }
            }
            if (!unique) {
                nduplicated++;
                constrains[i][0] = minConst - i;
                //printf("A duplicated at=%d\n",i);
            }
        }
        
    }

    /* Routine for tournament selection, it creates a new_pop from old_pop by performing tournament selection and the crossover */
    int tournamentSelection(int k) {
        //select first
        int sel = randomint(0,npop-1);
        int selected = sel;

        //printIndividual(genes[selected],ngene,objectives[selected],nobj,constrains[selected],nconst,0,0,selected, false);
        
        for (int i=1;i<k;i++) {
            sel = randomint(0,npop-1);
            
            if (this->bestTo(sel,this,selected)) {
                selected = sel;
                //printIndividual(genes[selected],ngene,objectives[selected],nobj,constrains[selected],nconst,0,0,selected, false);
            }
        }

        //printf("tournamentSelection DONE=%d\n",selected);
        
        return selected;
        
    }

    bool bestTo(int sel, Population *pop2, int selected) {
        if ((constrains[sel][0] < 0) && (pop2->constrains[selected][0] < 0)) {
            //both are infeasible, choose the less infeasible
            if (constrains[sel][0] > pop2->constrains[selected][0]) {
                return true;
            } else {
                return false;
            }
        } else if ((constrains[sel][0] < 0) && (pop2->constrains[selected][0] >= 0)) {
            //sel is infeasible, selected is not
            return false;
        } else if ((constrains[sel][0] >= 0) && (pop2->constrains[selected][0] < 0)) {
            //sel is feasible, selected is not, replace
            return true;
        } else  {
            //both are feasible, choose the best
            if (objectives[sel][0] > pop2->objectives[selected][0]) {
                return true;
            } else {
                return false;
            }
        }
    }


    int tournamentSelectionNSGA2(int k) {
        //select first
        int sel = randomint(0,npop-1);
        int selected = sel;
        
        for (int i=1;i<k;i++) {
            sel = randomint(0,npop-1);

            if (rank[sel] < rank[selected]) { //lower rank, better
                selected = sel;
            } else if (rank[sel] == rank[selected]) { //same rank, look for crowdDistance
                if (crowdDistance[sel] > crowdDistance[selected]) {
                    selected = sel;
                }
            }
        }
        //printf("Selected: %f\n",selected.fitness);
        
        return selected;
        
    }

    void gaussianMutation(T *individual,int ngene,float indpb,T minvalue,T maxvalue, float mu, float sigma) {
        T prev;
        T newval;
        
        for (int i=0;i<ngene;i++) {
            double rnd = randomf();
            if (rnd <= indpb) {
                prev = individual[i];
                newval = prev + randomf_gauss(mu, sigma);
                if (newval < minvalue) {
                    newval = minvalue;
                } else if (newval > maxvalue) {
                    newval = maxvalue;
                }

                //printf("Mutating gene[%d]-->[%d to %d]\n",i,prev,newval);
                individual[i] = newval;
            }
        }
    }

    void crossoverUniform(const T *parent1,const T *parent2,T *child1,T *child2,int ngene,float indpb) {
        for (int i=0;i<ngene;i++) {
            if (randomf() < indpb) {
                child1[i] = parent2[i];
                child2[i] = parent1[i];
            } else {
                child1[i] = parent1[i];
                child2[i] = parent2[i];
            }
        }
    }

    void calculateFrontRank(vector<vector<int>> &front) {
        int size = npop;
        
        //creating sets
        vector< vector< int > > S(size, vector< int >()); //solution dominated by i-th. create vector of size=size with empty sets
        vector<int> n(size,0); //domination count. create vector of size=size, initialized = 0

        front.clear();
        front.push_back(vector<int>());
        
        for (int p=0;p<size;p++) {
            for (int q=0;q<size;q++) {
                if (dominates(p,q)) {
                    S[p].push_back(q);
                } else if (dominates(q,p)) {
                    n[p]++;
                }
            }
            if (n[p] == 0) {
                rank[p] = 0;
                front[0].push_back(p);
                //printf("individual[%d] is rank 0\n",p);
            }
        }

        //CrowdDistance first front
        calculateCrowdingDistance(front[0]);
        
        int i=0;
        while (front[i].size() > 0) {
            
            vector<int> Q;
            for (int p : front[i]) {
                for (int q : S[p]) { //q has all dominated solutions by p
                    n[q]--;
                    if (n[q] == 0) {
                        rank[q] = i+1;
                        Q.push_back(q);
                    }
                }
            }
            //printf("size of next frontier[%d]=%d\n",i+1,Q.size());
            //CrowdDistance next front
            if (Q.size() > 0) {
                calculateCrowdingDistance(Q);
            }

            i++;
            front.push_back(Q);
        }
    }


    vector<vector<int>> selectFront(vector<vector<int>> &front,int newSize) {
        vector<vector<int>> newFront;
        int i = 0; //count of selection at
        int j = 0; //iterate over fronts
        while (i < newSize) {
            vector<int> &f = front[j];
            //printf("selecting front[%d]. front size=%d, current=%d\n",j,f.size(),i);
            if (f.size() <= newSize - i) {
                //this front must be fully copied
                i += f.size();
                newFront.push_back(vector<int>(f));
            } else {
                //this front must be partially copied
                newFront.push_back(vector<int>(f.begin(),f.begin() + (newSize - i)));
                break;
            }
            j++;
        }
        
        return newFront;
    }

    void calculateCrowdingDistance(vector<int> front) {
        int size = front.size();
        
        if (size == 1) {
            //printf("Frontier size 1: %d\n",front[0]);
            crowdDistance[front[0]] = inf;
        } else if (size == 2) {
            //printf("Frontier size 2: %d and %d\n",front[0],front[1]);
            crowdDistance[front[0]] = inf;
            crowdDistance[front[1]] = inf;
        } else if (size > 2) {
            double objectiveRange;
            //initialize all distance to 0
            for (int i=0;i<size;i++) {
                crowdDistance[front[i]] = 0.0;
            }

            for (int o=0;o<nobj;o++) {
                //fill objective o-th values
                vector< pair<double,int> > objectiveValue(size);
                for (int i=0;i<size;i++) {
                    objectiveValue[i] = make_pair(objectives[front[i]][o],front[i]);
                }
            
                //sort
                std::sort(objectiveValue.begin(), objectiveValue.end(),myfunction);
                objectiveRange = (objectiveValue[size-1].first - objectiveValue[0].first);

                crowdDistance[objectiveValue[0].second] = inf;
                crowdDistance[objectiveValue[size-1].second] = inf;
            
                for (int i=1;i<size-1;i++) {
                    if (objectiveRange != 0.0) {
                        double obj = (objectiveValue[i+1].first - objectiveValue[i-1].first)/objectiveRange;
                        //printf("OBJECTIVE[%d] Ind[%d][%d]. INC=%f. PREV=%f\n",o,i,objectiveValue[i].second,obj,pop->getPtr(objectiveValue[i].second)->crowdDistance);
                        crowdDistance[objectiveValue[i].second] += obj;
                    }
                }

                //for (int i=0;i<size;i++) {
                //    printf("OBJECTIVE[%d] Ind[%d][%d]=%f\n",o,i,front[i],pop->getPtr(front[i])->crowdDistance);
                //}

            }
            //for (int i=0;i<size;i++) {
            //    printf("FINALLY Ind[%d][%d]=%f\n",i,front[i],pop->getPtr(front[i])->crowdDistance);
            //}
        }
    }


    void SBXMutation(T *individual, int ngene, float indpb,float eta_m,T minvalue,T maxvalue) {
        int j;
        double delta,mut_pow;
        T y, newy;
        double range, rnd ;

        mut_pow = 1.0/(eta_m+1.0);
        range = double(maxvalue-minvalue);

        for (j=0; j<ngene; j++) {
            rnd = randomf();
            //printf("SBX %f < %f\n",rnd,indpb);
            if (rnd <= indpb) {
                y = individual[j];
                rnd = randomf();
                if (rnd <= 0.5)
                {
                    delta = pow(2.0*rnd,mut_pow) - 1.0;
                }
                else
                {
                    delta = 1.0 - pow(2.0*(1.0-rnd),mut_pow);
                }

                newy = y + delta*range;
                //printf("mutate %d delta=%f --> %f. min=%d, max=%d, rnd=%f\n",j,deltaq,y + deltaq*range,yl,yu,rnd);
                //printf("mutate[%d] %d -> %d. Min:Max=%d:%d\n",j,y,newy,minvalue,maxvalue);
                if (newy<minvalue)
                    newy = minvalue;
                if (newy>maxvalue)
                    newy = maxvalue;

                //printf("mutate[%d] %d -> %d. Min:Max=%d:%d\n",j,y,newy,minvalue,maxvalue);

                individual[j] = newy;
            }
        }
    }

    bool dominates(int a, int b) { 
        bool flag1 = false; //at least one objective is improved
        bool flag2 = false; //at least one objective is worse
        
        if (constrains[a][0] <0 && constrains[b][0] <0) { //both infeasible
            if (constrains[a][0] > constrains[b][0]) { //take the less infeasible
                return true;
            } else { //
                return false;
            }
        } else if (constrains[a][0] < 0 && constrains[b][0] >= 0) { //a is infeasible but b not
            return false;
        } else if (constrains[a][0] >= 0 && constrains[b][0] <0) { //a is feasible but b not
            return true;
        } else { //both are feasible
            for (int i=0; i<nobj; i++) {
                if (objectives[a][i] > objectives[b][i]) { //one objective is improved
                    flag1 = true;
                } else if (objectives[a][i] < objectives[b][i]) { //one objective is worse
                    flag2 = true;
                }
            } 
            if (flag1 && !flag2) {
                return true;
            //} else if (flag1==0 && flag2==1) {
            //    return false;
            } else {
                return false; //equals
            }
        }
    }

    void copyTo(int a, Population *other, int b) { 
        //for (int i=0;i<ngene;i++) {
        //    other->genes[b][i] = genes[a][i];
        //}
        memcpy(other->genes[b],genes[a],sizeof(T)*ngene);
        memcpy(other->objectives[b],objectives[a],sizeof(double)*nobj);
        memcpy(other->constrains[b],constrains[a],sizeof(double)*nconst);

        /*
        for (int i=0;i<nobj;i++) {
            other->objectives[b][i] = objectives[a][i];
        }
        for (int i=0;i<nconst;i++) {
            other->constrains[b][i] = constrains[a][i];
        }
        */
        
        other->rank[b] = rank[a];
        other->crowdDistance[b] = crowdDistance[a];
        other->evaluated[b] = evaluated[a];
    }


};


//using IntPopulation = Population<int>;

template class Population<int>;

#endif
