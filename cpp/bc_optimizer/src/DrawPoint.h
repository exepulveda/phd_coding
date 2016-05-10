#ifndef DRAWPOINT_H
#define DRAWPOINT_H

#include <vector>


#include "BlockModel.h"

#include <armadillo>

using namespace std;
using namespace arma;

class DrawPoint
{
    public:
        DrawPoint(BlockModel &bc);
        virtual ~DrawPoint();

        void setInfluence(Row<int> row);
        int extraction(int previousExtraction,int units, int *blocks);
        int extraction(int previousExtraction,int units, vector<int> &blocks);
        
        //int loc_i;
        //int loc_j;
        BlockModel &bm;

        int *stream;
        int streamsize;
        
        int size() { return streamsize; }
    protected:
    private:
        void computeInfluence();
};

#endif // DRAWPOINT_H
