#include <vector>
#include <set>
#include <stdio.h>

#include "Population.h"
#include "Individual.h"
#include "CrowdDistance.h"
#include "FrontRank.h"

using namespace std;

void calculateFrontRank(Population *population,vector<vector<int>> &front) {
    int size = population->size();
    
    //creating sets
    vector< vector< int > > S(size, vector< int >()); //solution dominated by i-th. create vector of size=size with empty sets
    vector<int> n(size,0); //domination count. create vector of size=size, initialized = 0

    Individual *indp;
    Individual *indq;

    front.clear();
    front.push_back(vector<int>());
    
    for (int p=0;p<size;p++) {
        indp = population->getPtr(p);
        for (int q=0;q<size;q++) {
            indq = population->getPtr(q);
            
            if (indp->dominates(indq)) {
                S[p].push_back(q);
            } else if (indq->dominates(indp)) {
                n[p]++;
            }
        }
        if (n[p] == 0) {
            indp->rank = 0;
            front[0].push_back(p);
            //printf("individual[%d] is rank 0\n",p);
        }
    }

    //CrowdDistance first front
    calculateCrowdingDistance(front[0],population);
    
    int i=0;
    while (front[i].size() > 0) {
        
        vector<int> Q;
        for (int p : front[i]) {
            for (int q : S[p]) { //q has all dominated solutions by p
                n[q]--;
                if (n[q] == 0) {
                    indq = population->getPtr(q);
                    indq->rank = i+1;
                    Q.push_back(q);
                }
            }
        }
        //printf("size of next frontier[%d]=%d\n",i+1,Q.size());
        //CrowdDistance next front
        if (Q.size() > 0) {
            calculateCrowdingDistance(Q,population);
        }

        i++;
        front.push_back(Q);
    }
}

void printAllFronts(Population *population,vector<vector<int>> &front) {
    int i=0;
    for (vector<int> &f : front) {
        printf("Front[%d], size=%d. Elements:\n",i++,f.size());
        printFront(population,f);
    }
}

void printFront(Population *population,vector<int> &f) {
    for (int p : f) {
        population->getPtr(p)->print(p);
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
