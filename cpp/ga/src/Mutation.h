#ifndef MUTATION_INCLUDE
#define MUTATION_INCLUDE

#include "Individual.h"
#include "Random.h"

void GaussianMutation(Individual &ind,float indpb,float minvalue,float maxvalue,float mu = 0.0, float sigma = 5.0);

template<class T>
void SBXMutation(Individual *individual, float indpb,float eta_m,T minvalue,T maxvalue) {
    int j;
    double delta,mut_pow;
    T y, newy;
    double range, rnd ;

    mut_pow = 1.0/(eta_m+1.0);
    range = double(maxvalue-minvalue);

    for (j=0; j<individual->size(); j++) {
        rnd = randomf();
        //printf("SBX %f < %f\n",rnd,indpb);
        if (rnd <= indpb) {
            y = individual->getInt(j);
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

            individual->set(j,newy);
        }
    }
}

#endif
