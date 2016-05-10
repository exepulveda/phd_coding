#include <algorithm>
#include <iostream>
#include <stdio.h>

#include "DrawPoint.h"


using namespace std;

DrawPoint::DrawPoint(BlockModel &_bm): bm(_bm) {
    this->stream = NULL;
}

DrawPoint::~DrawPoint() {
    if (stream != NULL) {
        delete [] stream;
    }
}


void DrawPoint::setInfluence(Row<int> blocks)
{
    int k;
    
    //find k unit -1
    for (k=0; k<blocks.n_cols;k++) {
        //printf("blockid[%d]=%d\n",k,blocks[k]);
        if (blocks[k] < 0) {
            break;
        }
    }
    this->stream = new int[k];
    for (int i=0; i<k;i++) {
        stream[i] = blocks[i];
    }
    
    streamsize = k;
}

int DrawPoint::extraction(int previousExtraction,int units, int *blocks)
{
    int extractionStart = previousExtraction;
    int extractionStop = extractionStart + units;

    if (extractionStop > size()) {
        extractionStop = size();
    }

    // cout << "extractionStart: " << extractionStart << endl;
    // cout << "extractionStop: " << extractionStop << endl;


    int n = extractionStop - extractionStart;
    
    if (n < 0) {
        printf("BAD THING, extractionStop=%d, extractionStart=%d, units=%d\n",extractionStop,extractionStart,units);
        exit(-99);
    }
    
    for (int i=extractionStart;i<extractionStop;i++) {
        blocks[i-extractionStart] = stream[i];
    }

    return n;

}

int DrawPoint::extraction(int previousExtraction,int units, vector<int> &blocks)
{
    int extractionStart = previousExtraction;
    int extractionStop = extractionStart + units;

    if (extractionStop > size()) {
        extractionStop = size();
    }

    // cout << "extractionStart: " << extractionStart << endl;
    // cout << "extractionStop: " << extractionStop << endl;


    int n = extractionStop - extractionStart;
    
    if (n < 0) {
        printf("BAD THING, extractionStop=%d, extractionStart=%d, units=%d\n",extractionStop,extractionStart,units);
        exit(-99);
    }
    
    blocks.clear();
    for (int i=extractionStart;i<extractionStop;i++) {
        blocks.push_back(stream[i]);
    }

    return n;

}
