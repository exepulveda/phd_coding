#include <stdio.h>
#include <assert.h>
#include <algorithm>
#include <cstring>

#include "Population.h"
#include "Random.h"
#include <boost/crc.hpp>

using namespace std;

template <typename T>
void Population<T>::initialize() {
    randomize();
            
    //first evaluation
    #pragma omp parallel for
    for (int i=0;i<npop;i++) {
        this->evaluator(genes[i],objectives[i],constrains[i]);
    }
}

template <typename T>
void Population<T>::evolve() {
    //one generation more
    //1.- create offspring
    //printf("creating offspring\n");
    //printf("(1) Best individual at the moment: %f\n",bestIndividual->fitness);
    //Population offspring = Population(this);
    //printf("(1.0) Best individual at the moment: %f\n",bestIndividual->fitness);
    
    Population<T> offspring = Population<T>(npop*2,ngene,nobj,nconst,evaluator);
    
    double rnd;
    for (int i=0;i<size();i+=2) {
        //select 2 parents by Tournoument selection
        T* parent1 = genes[select()];
        T* parent2 = genes[select()];
        //crossover parents
        T* child1 = offspring.genes[i];
        T* child2 = offspring.genes[i+1];
        
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
            memcpy(child1,parent1,sizeof(T)*ngene);
            memcpy(child2,parent2,sizeof(T)*ngene);
        }
        //printf("(1.1) Best individual at the moment: %f\n",bestIndividual->fitness);
    }
    //2.- mutate
    #pragma omp parallel for
    for (int i=0;i<npop;i++) {
        if (randomf() < mutationProb_) {
            //printf("mutation prob=%f\n",indProb_);
            //offspring.get(i).print();
            mutate(offspring.genes[i]);
            //offspring.get(i).print();
        }
    }
    //3.- Evaluate
    //printf("evaluation\n");
    //printf("(3) Best individual at the moment: %f\n",bestIndividual->fitness);
    #pragma omp parallel for
    for (int i=0;i<npop;i++) {
        //offspring.get(i).fitness = this->evaluator_(offspring.get(i));
        this->evaluator(offspring.genes[i],offspring.objectives[i],offspring.constrains[i]);
        //offspring.getRef(i).print();
    }

    
    //replace new generation
    //printf("replacing new generation\n");
    double mean_fit = 0;
    double min_fit = offspring.objectives[0][0];
    double max_fit = min_fit;
    
    for (int i=0;i<npop;i++) {
        memcpy(genes[i],offspring.genes[i],sizeof(T)*ngene);
        //if (bestIndividual->fitness[0] < individuals[i]->fitness[0]) {
        //    setNewBest(individuals[i]);
        //}
        mean_fit += objectives[i][0];
        if (min_fit > objectives[i][0]) {
            min_fit = objectives[i][0];
        }
        if (max_fit < objectives[i][0]) {
            max_fit = objectives[i][0];
        }
    }
    mean_fit /= npop;
    
    printf("FITNESS Min=%f | Mean=%f | Max=%f\n",min_fit,mean_fit,max_fit);
    
}

template <typename T>
int Population<T>::select() {
    return tournamentSelection(5);
}

template <typename T>
int Population<T>::selectNSGA2() {
    return tournamentSelectionNSGA2(5);
}


template <typename T>
void Population<T>::crossover(const T *parent1,const T *parent2,T *child1,T *child2) {
    crossoverUniform(parent1,parent2,child1,child2,ngene,0.5);
}

template <typename T>
void Population<T>::mutate(T *individual) {
    SBXMutation(individual,ngene,indProb_,10,minvalue_,maxvalue_);
}

template <typename T>
void Population<T>::purgeDuplicated() {
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
            printf("A duplicated at=%d\n",i);
        }
    }
    
}

/* Routine for tournament selection, it creates a new_pop from old_pop by performing tournament selection and the crossover */
template <typename T>
int Population<T>::tournamentSelection(int k) {
    //select first
    int sel = randomint(0,npop-1);
    int selected = sel;
    
    double max_fitness = objectives[sel][0];

    //printf("Selected: %d:%f\n",0,selected.fitness);

    
    for (int i=1;i<k;i++) {
        sel = randomint(0,npop-1);
        
        //printf("Select: %d:%f\n",i,ind.fitness);
        if (objectives[sel][0] > max_fitness) {
            max_fitness = objectives[sel][0];
            selected = sel;
        }
    }

    //printf("Selected: %f\n",selected.fitness);
    
    return selected;
    
}

template <typename T>
int Population<T>::tournamentSelectionNSGA2(int k) {
    //select first
    int sel = randomint(0,npop-1);
    int selected = sel;
    
    for (int i=1;i<k;i++) {
        sel = randomint(0,npop-1);

        if (rank[sel] < rank[selected]) { //lower rank, better
            selected = sel;
        } else if (rank[sel] == rank[selected]) { //same rank, look for creowDistance
            if (crowdDistance[sel] > crowdDistance[selected]) {
                selected = sel;
            }
        }
    }
    //printf("Selected: %f\n",selected.fitness);
    
    return selected;
    
}

template <typename T>
void Population<T>::gaussianMutation(T *individual,int ngene,float indpb,T minvalue,T maxvalue, float mu, float sigma) {
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

template <typename T>
void Population<T>::crossoverUniform(const T *parent1,const T *parent2,T *child1,T *child2,int ngene,float indpb) {
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

template <typename T>
void Population<T>::calculateFrontRank(vector<vector<int>> &front) {
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


template <typename T>
vector<vector<int>> Population<T>::selectFront(vector<vector<int>> &front,int newSize) {
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

bool myfunction(pair<double,int> i,pair<double,int> j) { return (i.first<j.first); }


template <typename T>
void Population<T>::calculateCrowdingDistance(vector<int> front) {
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


template <typename T>
void Population<T>::SBXMutation(T *individual, int ngene, float indpb,float eta_m,T minvalue,T maxvalue) {
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
            if (newy<minvalue)
                newy = minvalue;
            if (newy>maxvalue)
                newy = maxvalue;

            //printf("mutate %d %d -> %d\n",j,y,newy);

            individual[j] = newy;
        }
    }
}

template <typename T>
bool Population<T>::dominates(int a, int b) { 
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

template <typename T>
void Population<T>::copyTo(int a, Population *other, int b) { 
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
}