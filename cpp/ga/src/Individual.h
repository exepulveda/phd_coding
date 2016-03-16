#ifndef INDIVIDUAL_INCLUDE
#define INDIVIDUAL_INCLUDE

#include <cstddef>

class Individual;

class Individual {
    public:
        Individual(const size_t size) {
            fitness = 0;
            size_ = size;
        }

        virtual ~Individual() {};

    
        double fitness;
        
        virtual size_t size() { return size_; };
        
        virtual void set(int i, int value) = 0;
        virtual int getInt(int i) = 0;
        
        virtual void copy(Individual &other);
        virtual void copy(Individual *other);
        
        virtual void randomize(int minvalue,int maxvalue) = 0;
        
        virtual Individual *clone() = 0;
        virtual void print() = 0;

        
    protected:
        size_t size_;
};


#endif
