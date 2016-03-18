#include "Individual.h"

void Individual::copy(Individual &other) { 
    copy(&other);
}

void Individual::copy(Individual *other) { 
    fitness = other->fitness; 
    constrains_ = other->constrains_; 
    size_ = other->size_; 
}


bool Individual::dominates(Individual *b) { 
    int flag1; //at least one objective is improved
    int flag2; //at least one objective is worse
    flag1 = 0;
    flag2 = 0;
    Individual *a = this;
    
    if (a->constrains() <0 && b->constrains()<0) { //both infeasible
        if (a->constrains() > b->constrains()) { //take the less infeasible
            return true;
        } else { //
            return false;
        }
    } else if (a->constrains() < 0 && b->constrains() >= 0) { //a is infeasible but b not
        return false;
    } else if (a->constrains() >= 0 && b->constrains() <0) { //a is feasible but b not
        return true;
    } else { //both are feasible
        for (int i=0; i<nobj; i++) {
            if (a->fitness[i] > b->fitness[i]) {
                flag1 = 1;
            } else if (a->fitness[i] < b->fitness[i]) {
                flag2 = 1;
            }
        } if (flag1==1 && flag2==0) {
            return true;
        } else if (flag1==0 && flag2==1) {
            return false;
        } else {
            return false; //equals
        }
    }
}
