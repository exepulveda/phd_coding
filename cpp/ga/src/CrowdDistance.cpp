#include <algorithm>
#include "CrowdDistance.h"

/* Routine to compute crowding distance based on ojbective function values when the population in in the form of a list */

bool myfunction(pair<double,int> i,pair<double,int> j) { return (i.first<j.first); }


void calculateCrowdingDistance(vector<int> front,Population *pop) {
    int size = front.size();
    
    if (size == 1) {
        //printf("Frontier size 1: %d\n",front[0]);
        pop->getPtr(front[0])->crowdDistance = inf;
    } else if (size == 2) {
        //printf("Frontier size 2: %d and %d\n",front[0],front[1]);
        pop->getPtr(front[0])->crowdDistance = inf;
        pop->getPtr(front[1])->crowdDistance = inf;
    } else if (size > 2) {
        double objectiveRange;
        //initialize all distance to 0
        for (int i=0;i<size;i++) {
            pop->getPtr(front[i])->crowdDistance = 0.0;
        }

        int nobj = pop->getPtr(0)->nobj;

        for (int o=0;o<nobj;o++) {
            //fill objective o-th values
            vector< pair<double,int> > objectiveValue(size);
            for (int i=0;i<size;i++) {
                objectiveValue[i] = make_pair(pop->getPtr(front[i])->fitness[o],front[i]);
            }
        
            //sort
            sort(objectiveValue.begin(), objectiveValue.end(),myfunction);
            objectiveRange = (objectiveValue[size-1].first - objectiveValue[0].first);

            pop->getPtr(objectiveValue[0].second)->crowdDistance = inf;
            pop->getPtr(objectiveValue[size-1].second)->crowdDistance = inf;
        
            for (int i=1;i<size-1;i++) {
                if (objectiveRange != 0.0) {
                    double obj = (objectiveValue[i+1].first - objectiveValue[i-1].first)/objectiveRange;
                    printf("OBJECTIVE[%d] Ind[%d][%d]. INC=%f. PREV=%f\n",o,i,objectiveValue[i].second,obj,pop->getPtr(objectiveValue[i].second)->crowdDistance);
                    pop->getPtr(objectiveValue[i].second)->crowdDistance += obj;
                }
            }

            for (int i=0;i<size;i++) {
                printf("OBJECTIVE[%d] Ind[%d][%d]=%f\n",o,i,front[i],pop->getPtr(front[i])->crowdDistance);
            }

        }
        for (int i=0;i<size;i++) {
            printf("FINALLY Ind[%d][%d]=%f\n",i,front[i],pop->getPtr(front[i])->crowdDistance);
        }
    }
}
