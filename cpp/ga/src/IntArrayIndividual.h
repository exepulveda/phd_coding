#ifndef INTARRAYINDIVIDUAL_INCLUDE
#define INTARRAYINDIVIDUAL_INCLUDE

#include <array>
#include <cstddef>

#include "Individual.h"
#include "Random.h"

using namespace std;

class IntArrayIndividual: public Individual {
protected:

public:
    int *gene;
    IntArrayIndividual(const size_t N);

    virtual ~IntArrayIndividual();
      
    void set(int i, int value) {
        this->gene[i] = value;
    }

    int getInt(int i) {
        return gene[i];
    }
        
    void copy(Individual &other);
    void copy(Individual *other);

    Individual *clone() {
        IntArrayIndividual *newobj = new IntArrayIndividual(size());
        newobj->copy(*this);
        return newobj;
    }


    virtual void randomize(int minvalue,int maxvalue) {
        //printf("randomize - IntArrayIndividual\n");
        for (int i=0;i<size();i++) {
            gene[i] = randomint(minvalue,maxvalue);
            //printf("gene[%d]=%d\n",i,gene[i]);
        }       
    }
    
    void print() {
        int s=0;
        for (int i=0;i<size();i++) {
            printf("%d|",gene[i]);
            s += gene[i];
        }
        printf("fitness=%f|%d\n",fitness,s);
    }
    
};


#endif
