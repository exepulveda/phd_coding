#ifndef FRONTRANK_H
#define FRONTRANK_H

#include <vector>

#include "Population.h"

using namespace std;

void calculateFrontRank(Population *population,vector<vector<int>> &front);

void printAllFronts(Population *population,vector<vector<int>> &front);

void printFront(Population *population,vector<int> &f);

vector<vector<int>> selectFront(vector<vector<int>> &front,int size);

#endif
