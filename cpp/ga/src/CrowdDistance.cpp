#include <algorithm>
#include "CrowdDistance.h"

/* Routine to compute crowding distance based on ojbective function values when the population in in the form of a list */

void calculateCrowdingDistance(vector<int> front,Population *pop) {
    int size = front.size();
    
    //initialize all distance to 0
    vector<double> distances(size,0);
    
    vector<double> crowdDistance(size,0);
    vector< pair<double,int> > objectiveValue(size);
    double objectiveMax;
    double objectiveMin;
    double objectiveRange;
    
    int nobj = pop->getPtr(0)->nobj;

    for (int o=0;o<nobj;o++) {
        //fill objective o-th values
        if (size >= 2) {
            objectiveValue[0] = make_pair(pop->getPtr(front[0])->fitness[o],front[0]);
            objectiveMax = objectiveValue[0].first;
            objectiveMin = objectiveMax;
            for (int i=1;i<size;i++) {
                objectiveValue[i] = make_pair(pop->getPtr(front[i])->fitness[o],front[i]);
                if (objectiveMax < objectiveValue[i].first) {
                    objectiveMax = objectiveValue[i].first;
                }
                if (objectiveMin > objectiveValue[i].first) {
                    objectiveMin = objectiveValue[i].first;
                }
            }
            objectiveRange = objectiveMax - objectiveMin;
        
            //sort
            sort(objectiveValue.begin(), objectiveValue.end());
            if (size > 0) {
                crowdDistance[0] = inf;
                crowdDistance[size-1] = inf;
                
                pop->getPtr(0)->crowdDistance = crowdDistance[0];
                pop->getPtr(size-1)->crowdDistance = crowdDistance[size-1];
            }
        
            for (int i=1;i<size-1;i++) {
                crowdDistance[i] += (objectiveValue[i+1].first - objectiveValue[i-1].first)/objectiveRange;
                pop->getPtr(i)->crowdDistance = crowdDistance[i];
            }
        } else if (size == 1) {
            crowdDistance[0] = inf;
            pop->getPtr(0)->crowdDistance = crowdDistance[0];
        } else if (size == 2) {
            crowdDistance[0] = inf;
            crowdDistance[1] = inf;
            pop->getPtr(0)->crowdDistance = crowdDistance[0];
            pop->getPtr(1)->crowdDistance = crowdDistance[1];
        }
    }
    
}
