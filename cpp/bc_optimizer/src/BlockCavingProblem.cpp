#include <fstream>
#include <assert.h>
#include <set>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "BlockCavingProblem.h"
#include "DrawPoint.h"
#include "hdf5_utils.h"

namespace pt = boost::property_tree;

using namespace boost::property_tree;
using namespace pt::json_parser;

BlockCavingProblem::BlockCavingProblem()
{
    this->drawPoints = NULL;
}

BlockCavingProblem::~BlockCavingProblem()
{
    if (this->drawPoints) {
        for (int i=0;i<ndp;i++) {
            delete this->drawPoints[i];
        }

        delete [] this->drawPoints;
    }

    delete [] deduct;
    delete [] payableRate;
    delete [] refiningCharge;
    delete [] price;
    //delete [] grades;
    //delete [] concentrates;
    //delete [] recovery;
    //delete [] nsr;

}

void BlockCavingProblem::setupPeriods(int periods, float discountRate) {
    this->nperiods = periods;
    this->discount.resize(periods);

    for (int i=0;i<periods;i++) {
        this->discount(i) = 1.0/pow(1.0 + discountRate,i);
        printf("discount(%d)=%f\n",i,this->discount(i));
    }
}

void BlockCavingProblem::load(string filename)
{
    string datafile;
    string dataset;

    ptree pt;
    json_parser::read_json(filename, pt);


    int nx = pt.get<int>("blockModel.nodes.x");
    int ny = pt.get<int>("blockModel.nodes.y");
    int nz = pt.get<int>("blockModel.nodes.z");

    this->bm = BlockModel(nx,ny,nz);

    float sx = pt.get<float>("blockModel.sizes.x");
    float sy = pt.get<float>("blockModel.sizes.y");
    float sz = pt.get<float>("blockModel.sizes.z");

    this->bm.setVolume(sx,sy,sz);

    //periods
    //this->setupPeriods(root["periods"].asInt(),root["discount_rate"].asFloat());
    this->setupPeriods(pt.get<int>("periods"),pt.get<float>("discount_rate"));


    //drawpoints
    datafile = pt.get<string>("drawpoints.datafile"); //root["density"]["datafile"].asString();
    dataset = pt.get<string>("drawpoints.dataset"); // root["density"]["dataset"].asString();

    Mat<int> dpinfo;
    printf("loading drawpoints info...\n");
    load_hdf5_matrix< Mat< int > >(dpinfo, datafile, dataset);
    //printf("DP info:\n");
    //for (int i=0;i<5;i++) {
    //    for (int j=0;j<5;j++) {
    //        printf("DP[%d,%d]=%d\n",i,j,dpinfo(i,j));
    //   }
    //}    
    printf("DONE...\n");
    
    //imat dpi(dpinfo);
    
    
    //dpi = dpinfo;
    
    this->setupDrawPoints(dpinfo.n_rows, dpinfo);
    
    //dpi.clear();

    //load density to calculate tonnage
    datafile = pt.get<string>("density.datafile"); //root["density"]["datafile"].asString();
    dataset = pt.get<string>("density.dataset"); // root["density"]["dataset"].asString();

    load_hdf5_vector(tonnage, datafile, dataset);
    
    tonnage *= bm.getVolume();
    
    //for (int i=0;i<tonnage.n_elem;i++) {
    //   printf("tonnage[%d]=%f\n",i,tonnage[i]);
    //}

    //exit(0);

    //this->grades = new mat[n];
    //this->recovery = new mat[n];
    //this->concentrates = new mat[n];
    //MAT this->nsr_average = rowvec(bm.n);
    //MAT this->tonnage = rowvec(bm.n);

    //load nsr_average
    //datafile = root["nsr_average"]["datafile"].asString();
    //dataset = root["nsr_average"]["dataset"].asString();

    //load_hdf5_vector(nsr_average, datafile, dataset);

    //load nsr
    datafile = pt.get<string>("nsr.datafile"); //root["nsr"]["datafile"].asString();
    dataset = pt.get<string>("nsr.dataset"); //root["nsr"]["dataset"].asString();

    irowvec shape = hdf5_shape(datafile, dataset);
    
    if (shape.n_elem > 1) {
        load_hdf5_matrix(nsr, datafile, dataset);
        this->nsim = shape[1];
        printf("nsr:\n");
        for (int i=0;i<5;i++) {
            for (int j=0;j<5;j++) {
                printf("a[%d,%d]=%f\n",i,j,nsr(i,j));
            }
        }
        printf("checksum=%f\n",mean(mean(nsr)));
    } else {
        load_hdf5_vector(nsr_average, datafile, dataset);
        this->nsim = 0;
        printf("nsr:\n");
        for (int i=0;i<5;i++) {
            printf("a[%d]=%f\n",i,nsr_average(i));
        }
        printf("checksum=%f\n",mean(nsr_average));            
    }
    printf("nsim=%d\n",nsim);
    

    /*
    element = root["metals"]["cu"];

    string name = element["name"].asString();

    deduct[i] = element["deduct"].asFloat();
    payableRate[i] = element["payableRate"].asFloat();
    refiningCharge[i] = element["refiningCharge"].asFloat();
    refiningCostConvertion[i] = 1.0;
    price[i] = element["price"].asFloat();
    */

        //this->grades[i] = mat(bm.n,nsim);
        //this->concentrates[i] = mat(bm.n,nsim);
        //this->recovery[i] = mat(bm.n,nsim);
        //this->nsr[i] = mat(bm.n,nsim);


        //datafile = element["datafiles"]["grades"].asString();
        //grades[i].load(datafile,hdf5_binary);
        //inplace_trans(grades[i]);
        //assert(grades[i].n_rows == bm.n);
        //assert(grades[i].n_cols == nsim);

        //datafile = element["datafiles"]["concentrates"].asString();
        //concentrates[i].load(datafile,hdf5_binary);
        //inplace_trans(concentrates[i]);
        //assert(concentrates[i].n_rows == bm.n);
        //assert(concentrates[i].n_cols == nsim);

        //datafile = element["datafiles"]["nsr"].asString();
        //nsr[i].load(datafile,hdf5_binary);
        //inplace_trans(nsr[i]);
        //assert(nsr[i].n_rows == bm.n);
        //assert(nsr[i].n_cols == nsim);

        //recovery[i].fill(0.85);

    //}

    this->confidenceInterval = pt.get<float>("confidenceInterval"); //root["confidenceInterval"].asFloat();
    this->minTonnage = pt.get<float>("feed_production.minimum"); //root["feed_production"]["minimum"].asFloat();
    this->maxTonnage = pt.get<float>("feed_production.maximum"); //root["feed_production"]["maximum"].asFloat();
    this->targetProduction = pt.get<float>("target_production"); //root["target_production"].asFloat();
    this->units = pt.get<int>("units"); //root["units"].asInt();
    this->maxExtraction = pt.get<int>("max_block_extraction"); //root["units"].asInt();
    
    
    printf("Loading OK\n");
}
void BlockCavingProblem::setupDrawPoints(int ndp, Mat<int> & x)
{
    this->ndp = ndp;

    this->drawPoints = new DrawPoint*[ndp];

    maxsize = 0;

    for (int i=0;i<ndp;i++) {
        DrawPoint *dp = new DrawPoint(bm);

        dp->setInfluence(x.row(i));

        this->drawPoints[i] = dp;
        
        if (maxsize < dp->size()) {
            maxsize = dp->size();
        }
        
        //printf("DP[%d], influence size=%d\n",i,dp->size());
    }
    
    
}

DrawPoint *BlockCavingProblem::getDrawPoint(int index)
{
    assert (drawPoints != NULL);
    assert (index < ndp);
    return drawPoints[index];
}

double BlockCavingProblem::blockBenefit(double content) {
    return 0.0;
    // return (price - cost) * ( content * density); // * recovery);
}

double BlockCavingProblem::blockValue(double content) {
    return 0.0; //return price * content * density; // * recovery;
}

void BlockCavingProblem::npv(imat &schedule,double &objective)
{
    assert(schedule.n_rows == ndp);
    assert(schedule.n_cols == nperiods);

    DrawPoint *dp;

    int extractedBlocks;
    int nBlocksToExtract;
    double ton; 
    double npv_total,npv;
    unsigned int blockid;

    //check block is extracted once
    set<int> blockExtracted;

    int maxsize = 0;
    for (int i=0;i<ndp;i++) {
        dp = getDrawPoint(i);
        maxsize = maxsize < dp->size() ? dp->size(): maxsize;
    }
    
    int *blocks = new int[maxsize];

    for (int i=0;i<ndp;i++) {
        int prevExtractions = 0;

        dp = getDrawPoint(i);

        //cout << "processing DP: " << i<< endl;

        double totalTonnage = 0.0;

        for (int j=0;j<nperiods;j++) {

            nBlocksToExtract = units * schedule(i,j);
            extractedBlocks = dp->extraction(prevExtractions,nBlocksToExtract,blocks);

            //printf("dp[%d] at period[%d],prevExtractions=%d,nBlocksToExtract=%d,extractedBlocks=%d\n",i,j,prevExtractions,nBlocksToExtract,extractedBlocks);
            // cout << "processing period: " << j<< endl;

            for (int b=0;b<extractedBlocks;b++) {
                blockid = blocks[b];

                //printf("%d,%d,%d,%d\n",i,j,b,blockid);

                assert (blockExtracted.find(blockid) == blockExtracted.end());

                blockExtracted.insert(blockid);

                ton = tonnage[blockid];
                //cout << "processing extraction: " << b << endl;

                npv = this->nsr[blockid] * ton * discount[j] / 1000;
                npv_total += npv;

                //totalTonnage += ton;

                //productionPeriod(j) += ton;
            }
            prevExtractions += extractedBlocks;
        }
    }
    
    delete [] blocks;

    objective = npv_total;
}

void BlockCavingProblem::npv_cvar_npv(int *schedule,double objectives[],double constrains[]) {
    DrawPoint *dp;

    int extractedBlocks;
    int nBlocksToExtract;
    double npv_total,npv;
    unsigned int blockid;

    rowvec npvDistribution(nsim);
    npvDistribution.fill(0);
    
    rowvec ton(nperiods);
    ton.fill(0.0);
    
    //check block is extracted once
    int *blocks = new int[maxsize];
    int prevExtractions;
    double totalTonnageExtraction;
    double tmp_tonnage;
    double discounted_tonnage;
    //int index;
    //int i,j,b,k;
    

    //printf("maxsize=[%d]\n",maxsize);
    //printf("ndp=[%d]\n",ndp);
    //printf("nperiods=[%d]\n",nperiods);

    for (int i=0;i<ndp;i++) {
        prevExtractions = 0;
        dp = getDrawPoint(i);

        for (int j=0;j<nperiods;j++) {
            int index=i*nperiods +j;
            
            //if (index >= (nperiods * ndp)) {
            //    printf("Problem at %d,%d, index=%d. nperiods=%d,ndp=%d\n",i,j,index,nperiods,ndp);
            //    exit(0);
            //}
            
            nBlocksToExtract = units * schedule[index]; //[(i,j);
            extractedBlocks = dp->extraction(prevExtractions,nBlocksToExtract,blocks);
            schedule[index] = extractedBlocks;

            //printf("dp[%d] at period[%d],prevExtractions=%d,nBlocksToExtract=%d,extractedBlocks=%d\n",i,j,prevExtractions,nBlocksToExtract,extractedBlocks);
            // cout << "processing period: " << j<< endl;

            totalTonnageExtraction = 0.0;
            for (int b=0;b<extractedBlocks;b++) {
                blockid = blocks[b];
                tmp_tonnage = tonnage[blockid] / 1.0e3;
                totalTonnageExtraction += tmp_tonnage;

                //assert (blockExtracted.find(blockid) == blockExtracted.end());

                //blockExtracted.insert(blockid);
                discounted_tonnage = tmp_tonnage * discount[j] / 1.0e3;

                //npv_total += this->nsr_average[blockid] * discounted_tonnage;
                
                for (int k=0;k<nsim;k++) {
                    npvDistribution[k] += this->nsr(blockid,k) * discounted_tonnage;
                }
            }
            
            ton[j] += totalTonnageExtraction;
            //printf("%d,%d,%10.4f\n",i,j,totalTonnage);
            
            prevExtractions += extractedBlocks;
        }
    }

    //calculate mean NSR
    objectives[0] = mean(npvDistribution);
    //calculate cvar at risk level
    npvDistribution = sort(npvDistribution);
    //npvDistribution = cumsum(npvDistribution);
    int location_risk_level = int(nsim*(1.0-this->confidenceInterval)); //0.9
    
    objectives[1] = npvDistribution[location_risk_level] ;

    //constrains
    constrains[0] = 0.0;

    //computing if tonnage per period is in between min and max
    for (int j=0;j<nperiods;j++) {
        //printf("tonnage at period %d,min=%10.4f,max=%10.4f, tonnage= %10.4f\n",j,minTonnage,maxTonnage,ton[j]);
        if (ton[j] > maxTonnage) {
            constrains[0] += (maxTonnage - ton[j]);
        } else if (ton[j] < minTonnage) {
            constrains[0] += (ton[j] - minTonnage);
        }
    }

    delete [] blocks;
    
}

void BlockCavingProblem::calculateNSRTonnage(int *schedule,rowvec &npvDistribution,rowvec &tonnagePeriod) {
    double totalTonnageExtraction;
    double tmp_tonnage;
    double discounted_tonnage;
    int prevExtractions;
    int nBlocksToExtract;
    vector<int> blocks;
    DrawPoint *dp;

    npvDistribution.fill(0.0);
    tonnagePeriod.fill(0.0);

    //printf("calculateNSRTonnage ndp=%d,nperiods=%d,nsim=%d\n",ndp,nperiods,nsim);
    
    for (int i=0;i<ndp;i++) {
        
        dp = getDrawPoint(i);
        
        prevExtractions = 0;
        for (int p=0;p<nperiods;p++) {
            int index = i*nperiods + p;
            
            nBlocksToExtract = units * schedule[index];

            assert(nBlocksToExtract <= maxExtraction);
            assert(nBlocksToExtract >= 0);

            int efectiveExtraction = dp->extraction(prevExtractions,nBlocksToExtract,blocks);
            
            schedule[index] = efectiveExtraction; //fix

            
            totalTonnageExtraction = 0.0;
            for (int blockid : blocks) {
                //printf("dp[%d] at period[%d],prevExtractions=%d,nBlocksToExtract=%d,extractedBlocks=%d,blockid=%d\n",i,p,prevExtractions,nBlocksToExtract,blocks.size(),blockid);
                
                assert(blockid < bm.n);
                assert(blockid >= 0);

                tmp_tonnage = this->tonnage[blockid] / 1.0e3;
                totalTonnageExtraction += tmp_tonnage;
                discounted_tonnage = tmp_tonnage * discount[p] / 1.0e3;
                
                for (int k=0;k<nsim;k++) {
                    npvDistribution[k] += this->nsr(blockid,k) * discounted_tonnage;
                }
            }
            tonnagePeriod[p] += totalTonnageExtraction;
            
            prevExtractions += blocks.size();
        }
    }
    //printf("calculateNSRTonnage DONE\n");
}

void BlockCavingProblem::npv_cvar_npv(imat &schedule,double *objectives,double *constrains) {
    DrawPoint *dp;

    int extractedBlocks;
    int nBlocksToExtract;
    double npv_total,npv;
    unsigned int blockid;
    double discounted_tonnage;

    rowvec npvDistribution(nsim);
    npvDistribution.fill(0);
    
    rowvec ton(nperiods);
    ton.fill(0.0);
    
    //check block is extracted once
    int prevExtractions;
    double totalTonnage;
    double tmp_tonnage;

    //printf("maxsize=[%d]\n",maxsize);
    //printf("ndp=[%d]\n",ndp);
    //printf("nperiods=[%d]\n",nperiods);
    vector<int> blocks;

    for (int i=0;i<ndp;i++) {
        prevExtractions = 0;
        dp = getDrawPoint(i);
        
        for (int j=0;j<nperiods;j++) {

            nBlocksToExtract = units * schedule(i,j);
            extractedBlocks = dp->extraction(prevExtractions,nBlocksToExtract,blocks);
            schedule(i,j) = extractedBlocks;

            //printf("dp[%d] at period[%d],prevExtractions=%d,nBlocksToExtract=%d,extractedBlocks=%d\n",i,j,prevExtractions,nBlocksToExtract,extractedBlocks);
            // cout << "processing period: " << j<< endl;

            totalTonnage = 0.0;
            for (int b=0;b<extractedBlocks;b++) {
                blockid = blocks[b];
                tmp_tonnage = tonnage[blockid] / 1.0e3;
                totalTonnage += tmp_tonnage;

                //assert (blockExtracted.find(blockid) == blockExtracted.end());

                //blockExtracted.insert(blockid);
                discounted_tonnage = tmp_tonnage * discount[j] / 1.0e3;
                
                for (int k=0;k<nsim;k++) {
                    npv = this->nsr(blockid,k) * discounted_tonnage;
                    npvDistribution[k] += npv;
                }
            }
            
            ton[j] += totalTonnage;
            //printf("%d,%d,%10.4f\n",i,j,totalTonnage);
            
            prevExtractions += extractedBlocks;
        }
    }
    
    objectives[0] = mean(npvDistribution);
    //calculate cvar at risk level
    npvDistribution = sort(npvDistribution);
    //npvDistribution = cumsum(npvDistribution);
    int location_risk_level = int(nsim*(1.0-this->confidenceInterval)); //0.9
    
    
    objectives[1] = npvDistribution[location_risk_level];

    //constrains
    constrains[0] = 0.0;

    //computing if tonnage per period is in between min and max
    for (int j=0;j<nperiods;j++) {
        //printf("%d,min=%10.4f,max=%10.4f, tonnage= %10.4f\n",j,minTonnage,maxTonnage,ton(j));
        if (ton[j] > maxTonnage) {
            constrains[0] += (maxTonnage - ton[j]);
        } else if (ton[j] < minTonnage) {
            constrains[0] += (ton[j] - minTonnage);
        }
    }
    //if (constrains[0] >= 0.0) {
    //    printf("Feasible solution found\n");
    //}
   
}

void BlockCavingProblem::npv_tonage(imat &schedule,double *objectives)
{
    assert(schedule.n_rows == ndp);
    assert(schedule.n_cols == nperiods);

    DrawPoint *dp;

    int extractedBlocks;
    int nBlocksToExtract;
    double ton; 
    double npv_total,npv;
    unsigned int blockid;
    

    //check block is extracted once
    set<int> blockExtracted;

    int *blocks = new int[maxsize];

    double totalTonnage = 0.0;
    for (int i=0;i<ndp;i++) {
        int prevExtractions = 0;

        dp = getDrawPoint(i);

        //cout << "processing DP: " << i<< endl;


        for (int j=0;j<nperiods;j++) {

            nBlocksToExtract = units * schedule(i,j);
            extractedBlocks = dp->extraction(prevExtractions,nBlocksToExtract,blocks);

            //printf("dp[%d] at period[%d],prevExtractions=%d,nBlocksToExtract=%d,extractedBlocks=%d\n",i,j,prevExtractions,nBlocksToExtract,extractedBlocks);
            // cout << "processing period: " << j<< endl;

            for (int b=0;b<extractedBlocks;b++) {
                blockid = blocks[b];

                //printf("%d,%d,%d,%d\n",i,j,b,blockid);

                assert (blockExtracted.find(blockid) == blockExtracted.end());

                blockExtracted.insert(blockid);

                ton = tonnage[blockid];
                //cout << "processing extraction: " << b << endl;

                npv = this->nsr[blockid] * ton * discount[j] / 1000;
                npv_total += npv;

                totalTonnage += ton;

                //productionPeriod(j) += ton;
            }
            prevExtractions += extractedBlocks;
        }
    }
    
    delete [] blocks;

    objectives[0] = npv_total;
    objectives[1] = totalTonnage;
}

void BlockCavingProblem::single_npv_tonnage_deviation(int *schedule,double &nrs, double &deviation, int sim) {
    DrawPoint *dp;

    int extractedBlocks;
    int nBlocksToExtract;
    double tmp_tonnage,discounted_tonnage; 
    double npv_total,npv;
    unsigned int blockid;
    
    double *totalTonnage = new double[nperiods];
    int *blocks =  new int[maxsize];
    
    for (int p=0;p<nperiods;p++) {
        totalTonnage[p] = 0.0;
    }
    
    
    npv_total = 0.0;
    for (int i=0;i<ndp;i++) {
        int prevExtractions = 0;

        dp = getDrawPoint(i);

        //cout << "processing DP: " << i<< endl;


        for (int p=0;p<nperiods;p++) {
            int index = i*nperiods + p;

            nBlocksToExtract = units * schedule[index];
            extractedBlocks = dp->extraction(prevExtractions,nBlocksToExtract,blocks);
            
            schedule[index] = extractedBlocks; //fix

            //printf("dp[%d] at period[%d],prevExtractions=%d,nBlocksToExtract=%d,extractedBlocks=%d\n",i,p,prevExtractions,nBlocksToExtract,extractedBlocks);
            // cout << "processing period: " << j<< endl;

            for (int b=0;b<extractedBlocks;b++) {
                blockid = blocks[b];

                //assert(blockid < bm.n);
                //assert(blockid >= 0);

                tmp_tonnage = this->tonnage[blockid] / 1.0e3;
                totalTonnage[p] += tmp_tonnage;
                discounted_tonnage = tmp_tonnage * discount[p] / 1.0e3;

                if (sim < 0) {
                    npv = this->nsr_average[blockid] * discounted_tonnage;
                } else {
                    npv = this->nsr(blockid,sim) * discounted_tonnage;
                }
                npv_total += npv;

                //printf("dp[%d] at period[%d],blockid=%d,tonnage=%f,nsr_average=%f\n",i,p,blockid,this->tonnage[blockid],this->nsr_average[blockid]);

            }
            prevExtractions += extractedBlocks;
        }
    }
    
    nrs = npv_total;
    //deviation
    deviation = 0.0;
    for (int p=0;p<nperiods;p++) {
        //printf("deviation:: minTonnage=%f, maxTonnage=%f, totalTonnage[%d] =%f\n",minTonnage,maxTonnage,p+1,totalTonnage[p]);
        if (totalTonnage[p] < this->minTonnage) {
            deviation += (totalTonnage[p] - this->minTonnage);
        } else if (totalTonnage[p] > this->maxTonnage) {
            deviation += (this->maxTonnage - totalTonnage[p]);
        }
    }
    
    delete [] totalTonnage;
    delete [] blocks;

}

void BlockCavingProblem::average_npv_tonnage_deviation(int *schedule,double &nsr, double &deviation, double &constrain) {
    rowvec totalTonnage(this->nperiods);
    rowvec nsrSim(this->nsim);
    
    totalTonnage.fill(0.0);
    nsrSim.fill(0.0);
    
    calculateNSRTonnage(schedule,nsrSim,totalTonnage);
    
    nsr = mean(nsrSim);
    //deviation from production boundaries and target
    constrain = 0.0;
    deviation = 0.0;
    for (int p=0;p<nperiods;p++) {
        if (totalTonnage[p] < this->minTonnage) {
            constrain += (totalTonnage[p] - this->minTonnage);
        } else if (totalTonnage[p] > this->maxTonnage) {
            constrain += (this->maxTonnage - totalTonnage[p]);
        }
        //deviation
        deviation += fabs(totalTonnage[p] - this->targetProduction);
    }
    //deviation is the average over periods
    deviation = deviation / nperiods;
    

    /*
    if (deviation >= 0) {
        double m=0;
        for (int s=0;s<this->nsim;s++) {
            printf("%f|",nsr_sim(s));
            m += nsr_sim(s);
        }
        m = m / this->nsim;
        printf("|mean=%f\n",m);
    }
    */
}

void BlockCavingProblem::average_npv_variance(int *schedule,double &nsr, double &variance, double &constrain) {
    rowvec totalTonnage(this->nperiods);
    rowvec nsrSim(this->nsim);
    
    totalTonnage.fill(0.0);
    nsrSim.fill(0.0);
    
    calculateNSRTonnage(schedule,nsrSim,totalTonnage);
    
    //mean and variance at once
    nsr = mean(nsrSim);
    
    variance = 0.0;
    for (int i=0;i<nsim;i++) {
        variance += ((nsrSim[i] - nsr)*(nsrSim[i] - nsr)) / nsim;
    }

    //deviation from production boundaries and target
    constrain = 0.0;
    for (int p=0;p<nperiods;p++) {
        if (totalTonnage[p] < this->minTonnage) {
            constrain += (totalTonnage[p] - this->minTonnage);
        } else if (totalTonnage[p] > this->maxTonnage) {
            constrain += (this->maxTonnage - totalTonnage[p]);
        }
    }
}

void BlockCavingProblem::average_npv_cvar(int *schedule,double &nsr, double &cvar, double &constrain) {
    rowvec totalTonnage(this->nperiods);
    rowvec nsrSim(this->nsim);
    
    totalTonnage.fill(0.0);
    nsrSim.fill(0.0);
    
    calculateNSRTonnage(schedule,nsrSim,totalTonnage);
    
    //mean and variance at once
    nsr = mean(nsrSim);

    //cvar
    int location_risk_level = int(this->nsim*(1.0-this->confidenceInterval)); 
    nsrSim = sort(nsrSim);
    
    //cout << "nsrSim: " << nsrSim << endl;
    //cout << "location_risk_level: " << location_risk_level << endl;
    //cout << "nsr: " << nsr << ". cvar: " << cvar << endl;

    cvar = mean(nsrSim.subvec(0, location_risk_level));

    //deviation from production boundaries and target
    constrain = 0.0;
    for (int p=0;p<nperiods;p++) {
        if (totalTonnage[p] < this->minTonnage) {
            constrain += (totalTonnage[p] - this->minTonnage);
        } else if (totalTonnage[p] > this->maxTonnage) {
            constrain += (this->maxTonnage - totalTonnage[p]);
        }
    }
}


void BlockCavingProblem::average_nsr_var(int *schedule,double &nsr, double &nsrVar, double &constrain) {
    rowvec totalTonnage(this->nperiods);
    rowvec nsrSim(this->nsim);
    
    totalTonnage.fill(0.0);
    nsrSim.fill(0.0);
    
    calculateNSRTonnage(schedule,nsrSim,totalTonnage);
    
    //mean
    nsr = mean(nsrSim);
    //var
    nsrVar = var(nsrSim);
    
    //deviation from production boundaries and target
    constrain = 0.0;
    for (int p=0;p<nperiods;p++) {
        if (totalTonnage[p] < this->minTonnage) {
            constrain += (totalTonnage[p] - this->minTonnage);
        } else if (totalTonnage[p] > this->maxTonnage) {
            constrain += (this->maxTonnage - totalTonnage[p]);
        }
    }
    
}
/**
void BlockCavingProblem::npv_prod_deviation(imat &schedule,double *objectives)
{
    assert(schedule.n_rows == ndp);
    assert(schedule.n_cols == nperiods);

    DrawPoint *dp;

    urowvec blocks;

    double npv_average = 0;
    double npv_std;
    double grade_average = 0;
    double targetProduction = 6000.0;

    rowvec productionPeriod(nperiods);
    productionPeriod.fill(0.0);


    double bbs;
    double bb;
    double block_loss;
    int extractedBlocks;
    int nBlocksToExtract;
    double eb;

    rowvec ret(nsim);

    ret.fill(0.0);

    //check block is extracted once
    set<int> blockExtracted;

    for (int i=0;i<ndp;i++) {
        int prevExtractions = 0;

        dp = getDrawPoint(i);

        //cout << "processing DP: " << i<< endl;

        double totalTonnage = 0.0;

        for (int j=0;j<nperiods;j++) {

            nBlocksToExtract = units * schedule(i,j);
            extractedBlocks = dp->extraction(prevExtractions,nBlocksToExtract,blocks);

            // cout << "processing period: " << j<< endl;

            if (extractedBlocks > 0) {

                for (int b=0;b<extractedBlocks;b++) {
                    int blockid = blocks[b];

                    assert (blockExtracted.find(blockid) == blockExtracted.end());

                    blockExtracted.insert(blockid);

                    double ton = tonnage[blockid];
                    //cout << "processing extraction: " << b << endl;

                    for (int k=0;k<nsim;k++) {
                        //cout << "processing simulation: " << k << endl;


                        int nn = 2;

                        double npv=0.0;

                        for (int e=0;e<nn;e++) {
                            npv += this->nsr[e](blockid,k);
                        }

                        npv = npv * ton * discount(j) / 1000;

                        //printf("%d|%d|%d|%d|%f\n",i,j,blockid,k,npv);

                        ret(k) += npv;
                    }

                    totalTonnage += ton;

                    productionPeriod(j) += ton;
                }
                prevExtractions += extractedBlocks;
            }

        }
    }

    npv_average = mean(ret);
    npv_std = stddev(ret);

    //cout << "Statistics npv: " << npv_average << ":" << npv_std << endl;


    //exit(-0);
    //compute quantile
    ret = sort(ret);
    //exit(-0);
    //compute quantile
    //loss = sort(loss);
    double step = 1.0f/nsim;
    double s = 0;
    double w = 0;
    int i;


    ret = cumsum(ret);


    for (i=0;(i<nsim) && ((w+step) < 0.05);i++) {
        w += step;
        s = ret[i];
        //printf("step|w|min|max|s: %f|%f|%f|%f|%f\n",step,w,ret[0],ret[nsim-1],s);
    }

    //computing deviation from targetProduction
    double totalDeviation=0.0;
    double dev;
    for (int j=0;j<nperiods;j++) {
        dev = fabs(targetProduction - productionPeriod(j));
        totalDeviation += (dev / nperiods);
    }


    //s = i;

    //cout << "eval npv: " << npv_average << ":" << s << endl;

    //compute average
    //double average = mean(npv_ret);
    objectives[0] = npv_average;
    objectives[1] = totalDeviation;

}

void BlockCavingProblem::npv_distribution(imat &schedule,rowvec &npvDistribution,mat &tonnageDPPeriod)
{
    assert(schedule.n_rows == ndp);
    assert(schedule.n_cols == nperiods);

    DrawPoint *dp;
    urowvec blocks;

    rowvec ret(nsim);

    double npv_average = 0;
    double grade_average = 0;

    double bbs;
    double bb;
    double block_loss;
    int extractedBlocks;
    int nBlocksToExtract;
    double eb;

    npvDistribution.fill(0.0);
    tonnageDPPeriod.fill(0.0);


    for (int i=0;i<ndp;i++) {
        int prevExtractions = 0;

        dp = getDrawPoint(i);

        //cout << "processing DP: " << i<< endl;

        double totalTonnage = 0.0;

        for (int j=0;j<nperiods;j++) {

            nBlocksToExtract = units * schedule(i,j);
            extractedBlocks = dp->extraction(prevExtractions,nBlocksToExtract,blocks);

            double drawPointPeriodTonnage = 0.0;
            //cout << "processing period: " << j<< endl;

            if (extractedBlocks > 0) {

                for (int b=0;b<extractedBlocks;b++) {
                    int blockid = blocks[b];

                    double ton = tonnage[blockid];
                    //cout << "processing extraction: " << b << endl;

                    for (int k=0;k<nsim;k++) {
                        //cout << "processing simulation: " << k << endl;


                        int nn = 2;

                        double npv = 0.0;

                        for (int e=0;e<nn;e++) {
                            npv += this->nsr[e][blockid,k] * ton * discount(j) / 1000;
                        }

                        npv = npv * ton * discount(j) / 1000;

                        npvDistribution(k) += npv;

                        //printf("%d|%d|%d|%d|%f|%f|%f\n",i,j,blockid,k,npv,ton, discount(j));

                    }

                    totalTonnage += ton;

                    drawPointPeriodTonnage += ton;
                }
                prevExtractions += extractedBlocks;
            }

            tonnageDPPeriod(i,j) += drawPointPeriodTonnage;

        }
    }

    //exit(-0);
    //compute quantile
    npvDistribution = sort(npvDistribution);
}

void BlockCavingProblem::loss_distribution(imat &schedule,rowvec &ret)
{
    assert(schedule.n_rows == ndp);
    assert(schedule.n_cols == nperiods);

    DrawPoint *dp;

    vector<block_index> blocks;

    double npv_average = 0;
    double grade_average = 0;

    double bbs;
    double bb;
    double block_loss;
    int extractedBlocks;
    int nBlocksToExtract;
    double eb;

    ret.fill(0.0);


    for (int i=0;i<ndp;i++) {
        int prevExtractions = 0;

        dp = getDrawPoint(i);

        //cout << "processing DP: " << i<< endl;

        double totalTonnage = 0.0;

        for (int j=0;j<nperiods;j++) {

            blocks.clear();
            nBlocksToExtract = units * schedule(i,j);
            extractedBlocks = dp->extraction(prevExtractions,nBlocksToExtract,blocks);

            //cout << "processing period: " << j<< endl;

            if (extractedBlocks > 0) {

                for (int b=0;b<extractedBlocks;b++) {
                    int blockid = blocks[b];

                    double ton = tonnage[blockid];
                    //cout << "processing extraction: " << b << endl;

                    for (int k=0;k<nsim;k++) {
                        //cout << "processing simulation: " << k << endl;


                        int nn = 2;

                        double nvp=0.0;

                        for (int e=0;e<nn;e++) {
                            nvp += this->nsr[e][blockid,k] * ton * discount(j) / 1000;
                        }

                        //cout << i << ":" << j << ":" << b << ":" << nvp << endl;

                        ret(k) += nvp;
                    }

                    totalTonnage += ton;
                }
                prevExtractions += extractedBlocks;
            }

        }
    }

    //exit(-0);
    //compute quantile
    ret = sort(ret);
}
**/
