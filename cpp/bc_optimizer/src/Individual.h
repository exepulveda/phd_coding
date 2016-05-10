#ifndef INDIVIDUAL_INCLUDE
#define INDIVIDUAL_INCLUDE

#include <cstddef>
#include <vector>

using namespace std;

class Individual {
    public:
        Individual(const size_t size,int nobj=1,int nconst=0) {
            fitness = vector<double>(nobj,0);
            constrains_ = vector<double>(nconst,0);
            size_ = size;
            this->nconst = nconst;
            this->nobj = nobj;
        }

        virtual ~Individual() {
        }

    
        vector<double> fitness;
        vector<double> constrains_;
        
        virtual size_t size() { return size_; };
        
        virtual void set(int i, int value) = 0;
        virtual int getInt(int i) = 0;

        virtual int compare(Individual *other) = 0;
        
        virtual void copy(Individual &other);
        virtual void copy(Individual *other);
        
        virtual void randomize(int minvalue,int maxvalue) = 0;
        
        virtual Individual *clone() = 0;
        virtual void print(int tag = -1,bool printGene = false) = 0;

        virtual bool dominates(Individual *);

        virtual double constrains() { return constrains_[0]; }
        
        int size_;
        int nobj;
        int nconst;
        
        int rank;
        double crowdDistance;
};


#endif
