#ifndef BLOCKCAVINGPROBLEM_H
#define BLOCKCAVINGPROBLEM_H

#include <string>
#include <armadillo>

#include "BlockModel.h"
#include "DrawPoint.h"

using namespace arma;

using namespace std;

class BlockCavingProblem
{
    public:
        BlockCavingProblem();
        virtual ~BlockCavingProblem();

        void load(string filename);
        void setupDrawPoints(int ndp, Mat<int> & x);
        DrawPoint *getDrawPoint(int index);

        void npv(imat &schedule,double &objectives);
        void npv_tonage(imat &schedule,double *objectives);
        void npv_cvar_npv(imat &schedule,double *objectives,double *constrains);
        void npv_cvar_npv(int *schedule,double *objectives,double *constrains);
        
        //single objective
        void single_npv_tonnage_deviation(imat &schedule,double &nsr, double &deviation);
        void average_npv_tonnage_deviation(int *schedule,double &nsr, double &deviation, double &constrain);
        //multiobjective
        
        void average_npv_cvar(int *schedule,double &nsr, double &cvar, double &constrain);
        void average_nsr_var(int *schedule,double &nsr, double &nsrVar, double &constrain);
        void average_npv_variance(int *schedule,double &nsr, double &variance, double &constrain);
        
        void calculateNSRTonnage(int *schedule,rowvec &npvDistribution,rowvec &tonnage);
        
        //void npv_distribution(imat &schedule,rowvec &npvDistribution,mat &tonnage);
        //void npv_prod_deviation(imat &schedule,double *objectives);

        double blockBenefit(double content);
        double blockValue(double content);

        void setupPeriods(int periods, float discountRate);


        DrawPoint **drawPoints;
        int ndp;
        int nperiods;
        int nmetals;
        int units;
        int nsim;
        
        int maxExtraction;

        frowvec nsr_average;
        fmat nsr;
        frowvec tonnage;

        fmat fluorine_concentrate;
        
        int nblocks;

    protected:
    private:
        string filename;
        //mat grades;
        //mat recovery;
        //fmat concentrates;
        
        BlockModel bm;
        frowvec discount;

        float *deduct;
        float *payableRate;
        float *refiningCharge;
        float *refiningCostConvertion;
        float *price;

        float treatmentCharge;
        float confidenceInterval;
        
        float maxTonnage;
        float minTonnage;
        float targetProduction;
        
        int maxsize;

};

#endif // BLOCKCAVINGPROBLEM_H
