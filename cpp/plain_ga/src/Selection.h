#ifndef SELECTION_INCLUDE
#define SELECTION_INCLUDE


#include "Individual.h"
#include "Population.h"

int TournamentSelection(Population *population,int k);
int TournamentSelectionNSGA2(Population *population,int k);

#endif
