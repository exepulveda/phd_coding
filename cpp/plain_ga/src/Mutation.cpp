#include <stdio.h>

#include "Mutation.h"


void GaussianMutation(Individual &individual, float indpb,float minvalue,float maxvalue, float mu, float sigma) {
    int size = individual.size();
    int prev;
    int newval;
    
    for (int i=0;i<size;i++) {
        double rnd = randomf();
        if (rnd <= indpb) {
            prev = individual.getInt(i);
            newval = prev + randomf_gauss(mu, sigma);
            if (newval < minvalue) {
                newval = minvalue;
            } else if (newval > maxvalue) {
                newval = maxvalue;
            }

            printf("Mutating gene[%d]-->[%d to %d]\n",i,prev,newval);
            individual.set(i, newval);
        } else {
            //printf("No mutation prob %f < %f \n",rnd,indpb);
        }
    }
}

