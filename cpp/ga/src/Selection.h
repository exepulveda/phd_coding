#ifndef SELECTION_INCLUDE
#define SELECTION_INCLUDE


#include "Individual.h"
#include "Population.h"

Individual* TournamentSelection(Population *population,int k);
Individual* TournamentSelectionNSGA2(Population *population,int k);

#endif
