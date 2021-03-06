#ifndef INTARRAYINDIVIDUAL_INCLUDE
#define INTARRAYINDIVIDUAL_INCLUDE

#include <array>
#include <cstddef>
#include <cstring>

#include "Individual.h"
#include "Random.h"

using namespace std;

class IntArrayIndividual: public Individual {
protected:

public:
    int *gene;
    IntArrayIndividual(const size_t N,int nobj=1,int nconst=0);

    virtual ~IntArrayIndividual();
      
    void set(int i, int value) {
        this->gene[i] = value;
    }

    int getInt(int i) {
        return gene[i];
    }
        
    void copy(Individual &other);
    void copy(Individual *other);

    int compare(Individual *a) {
        IntArrayIndividual *aa = (IntArrayIndividual *)a;
        return memcmp(this->gene,aa->gene,sizeof(int)*size_);
    }


    Individual *clone() {
        IntArrayIndividual *newobj = new IntArrayIndividual(size(),nobj,nconst);
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
    
    void print(int tag=-1,bool printGene = false);
    
};


#endif
