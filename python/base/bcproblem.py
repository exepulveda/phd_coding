import h5py
import numpy as np
import json
import logging

logger = logging.getLogger(__name__)

from blockmodel import BlockModel
from drawpoint import DrawPoint

def cvar(distribution,risk_level=0.90):
    #print "deviation mean = ", deviation
    #quit()

    distribution.sort()
    n = distribution.shape[0]
    weights = np.empty(n)
    weights[:] = 1.0/n

    logger.debug("weight = %f",1.0/n)
    logger.debug("min deviation = %f",np.min(distribution))
    logger.debug("max deviation = %f",np.max(distribution))
    
    #find where is the percentil of risk measure, eg 95% confidence --> 0.95
    cummulate_probability = np.cumsum(weights)
    indices = np.where(cummulate_probability >= risk_level)[0]
    logger.debug("indices > 90 = %s",str(len(indices)))
    logger.debug("deviation > 90 = %s",str(distribution[indices]))
    #compute average 
    cvar = np.mean(distribution[indices])
    var = np.min(distribution[indices])
    logger.debug("cvar = %f",cvar)
    logger.debug("var = %f",var)

    return cvar,var
    

class BlockCavingProblem(object):
    
    def setup_periods(self,periods, discountRate):
        self.nperiods = periods
        self.discount = np.array([1.0/(1.0 + discountRate)**i for i in xrange(periods)])
        
    def setup_drawpoints(self,dp_blocks):
        self.ndp = len(dp_blocks)

        self.drawpoints = [None]*self.ndp

        for i in xrange(self.ndp):
            dp = DrawPoint(0, 0, self.bm)

            dp.set_influence(dp_blocks[i])
            #dp.set_influence(sizex, sizey, depth)

            self.drawpoints[i] = dp        

        #check for no duplicated blocks among DP
        for i in range(self.ndp):
            for j in range(i+1,self.ndp):
                duplicated  = set(self.drawpoints[i].stream) & set(self.drawpoints[j].stream)
                assert len(duplicated) == 0
        

    def load(self,filename,load_nsr=True):
        root = json.load(open(filename))

        node = root["blockModel"]

        nx,ny,nz = node["nodes"]["x"],node["nodes"]["y"],node["nodes"]["z"]
        sx,sy,sz = node["sizes"]["x"],node["sizes"]["y"],node["sizes"]["z"]
        
        nsim = int(root["nsim"])
        self.nsim = nsim
        self.bm = BlockModel(nx,ny,nz,sx,sy,sz);

        #metals
        #node = root["metals"]

        self.units = root["units"]
        self.max_block_extraction = root["max_block_extraction"]
        self.production_targets = root["target_production"]
        self.risk_level = root["confidenceInterval"]

        self.minimum_feed_production = root["feed_production"]["minimum"]
        self.maximum_feed_production = root["feed_production"]["maximum"]

        n = len(node)

        #dims = (nsim,self.bm.n)

        
        self.grades = None
        self.recovery = None
        self.concentrates = None
        self.nsr = None
        self.nsr_average = None
        
        
        self.deduct = np.empty(n)
        self.payableRate = np.empty(n)
        self.refiningCharge = np.empty(n)
        self.refiningCostConvertion = np.empty(n)
        self.price = np.empty(n)

        self.tonnage = np.empty(self.bm.n)

        nperiods = root["periods"]
        
        #loading Drawpoint
        try:
            dp_filename = root["drawpoints"]["datafile"]
            h5_dp = h5py.File(dp_filename, "r")
            dp_blocks_data = h5_dp[root["drawpoints"]["dataset"]][:,:]
            ndp,max_blocks = dp_blocks_data.shape
            h5_dp.close()
        except Exception as e:
            print "problem to open",dp_filename,root["drawpoints"]["dataset"]
            raise e
            
        dp_blocks = []
        for i in xrange(ndp):
            indices = np.where(dp_blocks_data[i,:] >= 0)[0]
            dp_blocks += [dp_blocks_data[i,indices]]

        self.setup_drawpoints(dp_blocks)
        self.setup_periods(nperiods,root["discount_rate"])        

        if load_nsr:
            nsr_filename = root["nsr"]["datafile"]
            h5_nsr = h5py.File(nsr_filename, "r")
            #load simulation
            self.nsr = np.array(h5_nsr[root["nsr"]["dataset"]])
            h5_nsr.close()

        #self.mining_cost = root["mining_cost"]

        element = root["density"]
        with h5py.File(element["datafile"],"r") as df:
            #density
            ds_name = element["dataset"]
            self.tonnage = df[ds_name][:] * self.bm.get_volume() # / 1000.0
            
        #logger.debug("dimension of production: %s",str(ndim))
        #self.production = self.production.flatten(axis=ndim[:-1])

        
        #logger.info("Total tonnage = %s",np.sum(self.tonnage))

        #self.reatmentCharge = float(root["treatmentCharge"])
        #self.nmetals = n

        #production main metal = Cu
        #self.production = self.tonnage * self.grades[0] * self.recovery[0] / self.concentrates[0] / 100.0
        #self.production = self.production.reshape(-1, self.production.shape[-1])
        #logger.debug("tonnage: %s",str(self.tonnage[0]))
        #logger.debug("grades: %s",str(self.grades[0,:,0]))
        #logger.debug("recovery: %s",str(self.recovery[0,:,:,0]))
        #logger.debug("concentrates: %s",str(self.concentrates[0,:,:,0]))
        #logger.debug("production shape: %s",str(self.production[:10,0]))
        
    def compute_objectives(self,schedule):
        assert(schedule.shape == (self.ndp,self.nperiods))

        npv_average = 0
        npv_std = 0
        grade_average = 0

        productionPeriod = np.zeros(self.nperiods)
        productionPeriod.fill(0.0);


        npv_sim = np.zeros(self.nsim)
        dev_sim = np.zeros(self.nsim)

        blockExtracted = set()

        for i in xrange(self.ndp):
            prevExtractions = 0

            dp = self.drawpoints[i]

            #//cout << "processing DP: " << i<< endl;

            totalTonnage = 0.0

            for j in xrange(self.nperiods):

                nBlocksToExtract = self.units * schedule[i,j]
                 
                blocks = dp.extraction(prevExtractions,nBlocksToExtract)
                extractedBlocks = len(blocks)
                prevExtractions += extractedBlocks

                #// cout << "processing period: " << j<< endl;

                #print i,j,extractedBlocks
                
                if (extractedBlocks > 0):

                    nsr_blocks = self.nsr[:,:,blocks]
                    ton_blocks = self.tonnage[blocks]
                    npv_sum = np.sum((nsr_blocks * ton_blocks) * (self.discount[j] ), axis=(0,2))
                    
                    #print nsr_blocks.shape
                    #print ton_blokcs.shape
                    #print npv_sum.shape

                    npv_sim += npv_sum
                    
                    
                    productionPeriod[j] += np.sum(ton_blocks)


        dev = np.abs(productionPeriod - self.targetProduction)

        return np.mean(npv_sim),np.sum(dev)

    def calculate_extracted_blocks(self,schedule,opening,closing):
        n = 0
        for i in xrange(self.ndp):
            prevExtractions = 0
            dp = self.drawpoints[i]
            opening_period = opening[i]
            closing_period = closing[i]

            for j in xrange(opening_period,closing_period+1):
                nBlocksToExtract = self.units * schedule[i,j]
                blocks = dp.extraction(prevExtractions,nBlocksToExtract)
                extractedBlocks = len(blocks)
                schedule[i,j] = extractedBlocks
                prevExtractions += extractedBlocks
                
                n +=  extractedBlocks

        return n

    def calculate_production_concentartes_distribution(self,schedule,opening,closing):
        assert(schedule.shape == (self.ndp,self.nperiods))
        assert(schedule.shape[0] == opening.shape[0])
        assert(schedule.shape[0] == closing.shape[0])

        productionPeriod = np.zeros(self.nperiods)

        production_per_period = np.zeros((self.production.shape[0],self.nperiods))

        for i in xrange(self.ndp):
            prevExtractions = 0

            dp = self.drawpoints[i]

            logger.debug("processing DP: %d",i)

            totalTonnage = 0.0


            opening_period = opening[i]
            closing_period = closing[i]

            #print i,opening_period,closing_period

            for j in xrange(opening_period,closing_period+1):
                nBlocksToExtract = self.units * schedule[i,j]
                 
                blocks = dp.extraction(prevExtractions,nBlocksToExtract)
                extractedBlocks = len(blocks)
                prevExtractions += extractedBlocks

                if (extractedBlocks > 0):
                    '''
                    production = tonnage * dmt
                    dmt = ore_grade * recovery / concentrate_grade / 100.0
                    
                    production can be calculated as input
                    '''                   
                    production = np.sum(self.production[:,blocks],axis=1)

                    #print i,j,schedule[i,j],production.shape,np.mean(production)

                    #logger.debug("production shape for (%d,%d) = %s",i,j,production.shape)

                    
                    production_per_period[:,j] += production

        
        return production_per_period

    def calculate_production_deviation_cvar(self,schedule,opening=None,closing=None):
        assert(schedule.shape == (self.ndp,self.nperiods))

        if opening is None:
            opening = np.zeros(self.ndp,dtype=np.int)
        else:
            assert(schedule.shape[0] == opening.shape[0])
            
        if closing is None:
            closing = np.zeros(self.ndp,dtype=np.int)
            closing[:] = self.nperiods - 1
        else:
            assert(schedule.shape[0] == closing.shape[0])

        productionPeriod = np.zeros(self.nperiods)

        production_per_period = np.zeros((self.production.shape[0],self.nperiods))

        for i in xrange(self.ndp):
            prevExtractions = 0

            dp = self.drawpoints[i]

            logger.debug("processing DP: %d",i)

            opening_period = opening[i]
            closing_period = closing[i]

            #print i,opening_period,closing_period

            for j in xrange(opening_period,closing_period+1):
                nBlocksToExtract = self.units * schedule[i,j]
                 
                blocks = dp.extraction(prevExtractions,nBlocksToExtract)
                extractedBlocks = len(blocks)
                prevExtractions += extractedBlocks

                if (extractedBlocks > 0):
                    '''
                    production = tonnage * dmt
                    dmt = ore_grade * recovery / concentrate_grade / 100.0
                    
                    production can be calculated as input
                    '''                   
                    production = np.sum(self.production[:,blocks],axis=1)
                    
                    #logger.debug("production shape for (%d,%d) = %s",i,j,production.shape)

                    
                    production_per_period[:,j] += production

        #print "production mean = ", np.mean(production_per_period,axis=0)

        #np.savetxt("production.csv",production_per_period)
        #quit()
        #deviation has the distribution of deviations in all periods
        #logger.info("deviation mean = %s",str(np.mean(production_per_period)))
        deviation = np.sum(np.abs(production_per_period - self.production_targets),axis=1)
        #print "deviation mean = ", deviation
        #quit()

        logger.debug("deviation shape = %s",str(deviation.shape))
        deviation.sort()
        nsim = production_per_period.shape[0]
        weights = np.empty(nsim)
        weights[:] = 1.0/nsim

        logger.debug("weight = %f",1.0/nsim)
        logger.debug("min deviation = %f",np.min(deviation))
        logger.debug("max deviation = %f",np.max(deviation))
        
        #find where is the percentil of risk measure, eg 95% confidence --> 0.95
        cummulate_probability = np.cumsum(weights)
        indices = np.where(cummulate_probability >= self.risk_level)[0]
        logger.debug("indices > 90 = %s",str(len(indices)))
        logger.debug("deviation > 90 = %s",str(deviation[indices]))
        #compute average 
        var = np.min(deviation[indices])
        cvar = np.mean(deviation[indices])
        logger.debug("cvar = %f",cvar)

        return cvar,var

    def calculate_npv_production(self,schedule,opening,closing):
        assert(schedule.shape == (self.ndp,self.nperiods))
        assert(schedule.shape[0] == opening.shape[0])
        assert(schedule.shape[0] == closing.shape[0])

        productionPeriod = np.zeros(self.nperiods)

        production_per_period = np.zeros((self.production.shape[0],self.nperiods))
        
        npv_sim = 0.0
        
        for i in xrange(self.ndp):
            prevExtractions = 0

            dp = self.drawpoints[i]

            #logger.debug("processing DP: %d",i)

            totalTonnage = 0.0


            opening_period = opening[i]
            closing_period = closing[i]

            #print i,opening_period,closing_period

            for j in xrange(opening_period,closing_period+1):
                nBlocksToExtract = self.units * schedule[i,j]
                 
                blocks = dp.extraction(prevExtractions,nBlocksToExtract)
                extractedBlocks = len(blocks)
                prevExtractions += extractedBlocks

                if (extractedBlocks > 0):
                    '''
                    production = tonnage * dmt
                    dmt = ore_grade * recovery / concentrate_grade / 100.0
                    
                    production can be calculated as input
                    '''                   
                    production = np.sum(self.production[:,blocks],axis=1)
                    
                    #logger.debug("production shape for (%d,%d) = %s",i,j,production.shape)

                    
                    production_per_period[:,j] += production
                    
                    '''npv'''
                    nsr_blocks = self.nsr_average[blocks]
                    ton_blocks = self.tonnage[blocks]
                    npv_sum = np.sum((nsr_blocks * ton_blocks) * (self.discount[j] ))
                    
                    #print nsr_blocks.shape
                    #print ton_blokcs.shape
                    #print npv_sum.shape

                    npv_sim += npv_sum


        #print "production mean = ", np.mean(production_per_period,axis=0)

        #np.savetxt("production.csv",production_per_period)
        #quit()
        #deviation has the distribution of deviations in all periods
        deviation = np.sum(np.abs(production_per_period - self.production_targets),axis=1)
        #print "deviation mean = ", deviation
        #quit()

        logger.debug("deviation shape = %s",str(deviation.shape))
        deviation.sort()
        nsim = production_per_period.shape[0]
        weights = np.empty(nsim)
        weights[:] = 1.0/nsim

        logger.debug("weight = %f",1.0/nsim)
        logger.debug("min deviation = %f",np.min(deviation))
        logger.debug("max deviation = %f",np.max(deviation))
        
        #find where is the percentil of risk measure, eg 95% confidence --> 0.95
        cummulate_probability = np.cumsum(weights)
        indices = np.where(cummulate_probability >= self.risk_level)[0]
        logger.debug("indices > 90 = %s",str(len(indices)))
        logger.debug("deviation > 90 = %s",str(deviation[indices]))
        #compute average 
        cvar = np.mean(deviation[indices])
        logger.debug("cvar = %f",cvar)

        return (npv_sim,cvar)

    def calculate_maximum_nsr(self,schedule,nsr=None):
        if nsr is None:
            nsr_data = self.nsr_average
        else:
            nsr_data = nsr

        assert (len(nsr_data.shape) == 1)

        npv = 0.0
        
        tonnage = np.zeros(self.nperiods)

        for i in xrange(self.ndp):
            prevExtractions = 0

            dp = self.drawpoints[i]
            for j in xrange(0,self.nperiods):
                nBlocksToExtract = self.units * schedule[i,j]
                 
                blocks = dp.extraction(prevExtractions,nBlocksToExtract)
                extractedBlocks = len(blocks)
                prevExtractions += extractedBlocks
                schedule[i,j] = extractedBlocks

                if (extractedBlocks > 0):
                    '''npv'''
                    nsr_blocks = nsr_data[blocks]
                    ton_blocks = self.tonnage[blocks] / 1.0e3
                    npv_sum = np.sum((nsr_blocks * ton_blocks)) * (self.discount[j] / 1000.0)
                    #npv_sum = np.sum(((nsr_blocks - self.mining_cost) * ton_blocks))
                    
                    npv += npv_sum
                    
                    tonnage[j] += np.sum(ton_blocks)


        #print np.min(tonnage),np.max(tonnage),self.minimum_feed_production,self.maximum_feed_production
        indices_less = np.where(tonnage < self.minimum_feed_production)[0]
        indices_more = np.where(tonnage > self.maximum_feed_production)[0]

        if len(indices_less) > 0 or len(indices_more) > 0:
            npv = np.sum(tonnage[indices_less] - self.minimum_feed_production) - np.sum(tonnage[indices_more] - self.maximum_feed_production)

        return npv,

    def calculate_maximum_average_nsr(self,schedule,nsr=None):
        if nsr is None:
            nsr_data = self.nsr
        else:
            nsr_data = nsr

        #nsr values must be 2D
        assert (len(nsr_data.shape) == 2)
        n,nsim = nsr_data.shape

        nsr = np.zeros(nsim)
        
        tonnage = np.zeros(self.nperiods)

        for i in xrange(self.ndp):
            prevExtractions = 0

            dp = self.drawpoints[i]
            for j in xrange(0,self.nperiods):
                nBlocksToExtract = self.units * schedule[i,j]
                 
                blocks = dp.extraction(prevExtractions,nBlocksToExtract)
                extractedBlocks = len(blocks)
                prevExtractions += extractedBlocks
                schedule[i,j] = extractedBlocks

                if (extractedBlocks > 0):
                    ton_blocks = self.tonnage[blocks] / 1.0e3
                    tonnage[j] += np.sum(ton_blocks)

                    '''npv'''
                    nsr_blocks = nsr_data[blocks,:]
                    nsr_sum = np.sum((nsr_blocks.T * ton_blocks),axis=1) * (self.discount[j] / 1.0e3)
                    #npv_sum = np.sum(((nsr_blocks - self.mining_cost) * ton_blocks))
                        
                    nsr += nsr_sum
                    
        #print np.min(tonnage),np.max(tonnage),self.minimum_feed_production,self.maximum_feed_production
        indices_less = np.where(tonnage < self.minimum_feed_production)[0]
        indices_more = np.where(tonnage > self.maximum_feed_production)[0]

        if len(indices_less) > 0 or len(indices_more) > 0:
            tonnage_deviation = np.sum(tonnage[indices_less] - self.minimum_feed_production) - np.sum(tonnage[indices_more] - self.maximum_feed_production)

            return tonnage_deviation,
        else:
            return np.mean(nsr),

    def calculate_npv_tonnage(self,schedule,nsr=None):
        if nsr is None:
            nsr_data = self.nsr_average
        else:
            nsr_data = nsr

        tonnage = np.zeros(self.nperiods)
        npv = np.zeros(self.nperiods)

        for i in xrange(self.ndp):
            prevExtractions = 0

            dp = self.drawpoints[i]
            totalTonnage = 0.0
            for j in xrange(0,self.nperiods):
                nBlocksToExtract = self.units * schedule[i,j]
                 
                blocks = dp.extraction(prevExtractions,nBlocksToExtract)
                extractedBlocks = len(blocks)
                prevExtractions += extractedBlocks
                schedule[i,j] = extractedBlocks

                if (extractedBlocks > 0):
                    '''npv'''
                    nsr_blocks = nsr_data[blocks]
                    ton_blocks = self.tonnage[blocks] / 1.0e3
                    npv_sum = np.sum((nsr_blocks * ton_blocks)) * (self.discount[j] / 1000.0)
                    #npv_sum = np.sum(((nsr_blocks - self.mining_cost) * ton_blocks))
                    
                    npv[j] += npv_sum
                    
                    tonnage[j] += np.sum(ton_blocks)


        return npv,tonnage

    def calculate_average_npv_tonnage(self,schedule,nsr=None):
        if nsr is None:
            nsr_data = self.nsr_average
        else:
            nsr_data = nsr

        tonnage = np.zeros(self.nperiods)
        npv = np.zeros(self.nperiods)
        
        n,nsim = nsr.shape
        
        for i in xrange(self.ndp):
            prevExtractions = 0

            dp = self.drawpoints[i]
            totalTonnage = 0.0
            for j in xrange(0,self.nperiods):
                nBlocksToExtract = self.units * schedule[i,j]
                 
                blocks = dp.extraction(prevExtractions,nBlocksToExtract)
                extractedBlocks = len(blocks)
                prevExtractions += extractedBlocks
                schedule[i,j] = extractedBlocks

                if (extractedBlocks > 0):
                    '''npv'''
                    nsr_blocks = np.mean(nsr_data[blocks],axis=1)
                    ton_blocks = self.tonnage[blocks] / 1.0e3
                    npv_sum = np.sum((nsr_blocks * ton_blocks)) * (self.discount[j] / 1000.0)
                    #npv_sum = np.sum(((nsr_blocks - self.mining_cost) * ton_blocks))
                    
                    npv[j] += npv_sum
                    
                    tonnage[j] += np.sum(ton_blocks)


        return npv,tonnage

    def calculate_all_npv_tonnage(self,schedule,nsr=None):
        if nsr is None:
            nsr_data = self.nsr_average
        else:
            nsr_data = nsr

        n,nsim = nsr.shape

        tonnage = np.zeros(self.nperiods)
        npv = np.zeros((nsim,self.nperiods))
        
        
        for i in xrange(self.ndp):
            prevExtractions = 0

            dp = self.drawpoints[i]
            totalTonnage = 0.0
            for j in xrange(0,self.nperiods):
                nBlocksToExtract = self.units * schedule[i,j]
                 
                blocks = dp.extraction(prevExtractions,nBlocksToExtract)
                extractedBlocks = len(blocks)
                prevExtractions += extractedBlocks
                schedule[i,j] = extractedBlocks

                if (extractedBlocks > 0):
                    '''npv'''
                    nsr_blocks = nsr_data[blocks].T
                    ton_blocks = self.tonnage[blocks] / 1.0e3
                    npv_sum = np.sum((ton_blocks * nsr_blocks),axis=1) * (self.discount[j] / 1000.0)
                    #npv_sum = np.sum(((nsr_blocks - self.mining_cost) * ton_blocks))
                    
                    npv[:,j] += npv_sum
                    
                    tonnage[j] += np.sum(ton_blocks)


        return npv,tonnage

    def calculate_plain_npv(self):
        '''npv'''
        nsr_blocks = self.nsr_average[:]
        ton_blocks = self.tonnage[:]
        #npv_sum = np.sum((nsr_blocks * ton_blocks) * (self.discount[j] / 1000))
        npv_sum = np.sum((nsr_blocks * ton_blocks))
        
        return npv_sum


    def calculate_npv_dp(self,i,schedule,opening_period,closing_period):
        npv_sim = 0.0
        
        prevExtractions = 0

        dp = self.drawpoints[i]

        #logger.debug("processing DP: %d",i)

        totalTonnage = 0.0

        for j in xrange(opening_period,closing_period+1):
            nBlocksToExtract = self.units * schedule[i,j]
             
            blocks = dp.extraction(prevExtractions,nBlocksToExtract)
            extractedBlocks = len(blocks)
            prevExtractions += extractedBlocks

            if (extractedBlocks > 0):
                '''npv'''
                nsr_blocks = self.nsr_average[blocks]
                ton_blocks = self.tonnage[blocks]
                npv_sum = np.sum((nsr_blocks* ton_blocks) * (self.discount[j]))
                
                #print nsr_blocks.shape
                #print ton_blokcs.shape
                #print npv_sum.shape


                npv_sim += npv_sum

        return npv_sim


    def calculate_dp_stats(self,i):
        dp = self.drawpoints[i]

        blocks = dp.stream[:]
        extractedBlocks = len(blocks)

        '''npv'''
        nsr_blocks = self.nsr_average[blocks]
        ton_blocks = self.tonnage[blocks]
        nsr_sum = np.sum(nsr_blocks * ton_blocks)
        
        #production average over simulations
        prod_mean = np.mean(self.production[:,blocks],axis=0)
        prod_min = np.min(self.production[:,blocks],axis=0)
        prod_max = np.max(self.production[:,blocks],axis=0)
        
        prod_dev = np.std(np.sum(self.production[:,blocks],axis=0))
                
        return nsr_sum,np.sum(ton_blocks),np.sum(prod_mean),np.sum(prod_min),np.sum(prod_max),prod_dev

    def calculate_npv_constrained(self,schedule,opening=None,closing=None):
        assert(schedule.shape == (self.ndp,self.nperiods))
        
        if opening is None:
            opening = np.zeros(self.ndp,dtype=np.int)
        else:
            assert(schedule.shape[0] == opening.shape[0])
            
        if closing is None:
            closing = np.zeros(self.ndp,dtype=np.int)
            closing[:] = self.nperiods - 1
        else:
            assert(schedule.shape[0] == closing.shape[0])

        npv = np.zeros(self.ndp)
        
        for i in xrange(self.ndp):
            npv[i] = self.calculate_npv_dp(i,schedule,opening[i],closing[i])

        return npv

    def calculate_npv_constrained_proxy(self,proxy,schedule,opening,closing):
        i,opening_index,closing_index,j,k,extraction = proxy.perturbation
        '''only DP i and j change, therefore only NVP calculation of these DP is needed
        '''

        schedule_tmp = schedule.copy()
        schedule_tmp[j,k] = extraction
        
        opening_tmp = opening.copy()
        closing_tmp = closing.copy()
        
        opening_tmp[i] = opening_index
        closing_tmp[i] = closing_index
        
        pnpv = np.sum(proxy._solution.npv)
        
        proxy._solution.npv[i] = self.calculate_npv_dp(i,schedule_tmp,opening_tmp[i],opening_tmp[i])
        proxy._solution.npv[j] = self.calculate_npv_dp(j,schedule_tmp,opening_tmp[j],opening_tmp[j])

        #print "previous npv",pnpv,"new npv",np.sum(proxy._solution.npv)

        
        return proxy._solution.npv

    def calculate_feed_tonnage_constrain(self,schedule,opening = None,closing = None):
        if opening is None:
            opening = np.zeros(self.ndp,dtype=np.int)
        else:
            assert(schedule.shape[0] == opening.shape[0])
            
        if closing is None:
            closing = np.zeros(self.ndp,dtype=np.int)
            closing[:] = self.nperiods - 1
        else:
            assert(schedule.shape[0] == closing.shape[0])


        production_period = self.calculate_feed_tonnage(schedule,opening,closing)

        #calculate the deviation from feed targets
        #logger.debug("minimum_feed_production=%f",self.minimum_feed_production)
        #logger.debug("maximum_feed_production=%f",self.maximum_feed_production)

        minp = np.zeros_like(production_period)
        indices = np.where(production_period < self.minimum_feed_production)[0]
        if len(indices) > 0:
            minp[indices] = self.minimum_feed_production - production_period[indices]

        maxp = np.zeros_like(production_period)
        indices = np.where(production_period > self.maximum_feed_production)[0]
        if len(indices) > 0:
            maxp[indices] = production_period[indices] - self.maximum_feed_production

            
        return tuple(maxp) + tuple(minp)

    def calculate_feed_tonnage(self,schedule,opening=None,closing=None):
        production_period = np.zeros(self.nperiods)

        for i in xrange(self.ndp):
            prevExtractions = 0

            dp = self.drawpoints[i]

            #logger.debug("processing DP: %d",i)

            if opening is None:
                opening_period=0
            else:
                opening_period = opening[i]
                
            if closing is None:
                closing_period=self.nperiods - 1
            else:
                closing_period = closing[i]

            #print i,opening_period,closing_period
            #logger.debug("DP[%d]: %d/%d",i+1,opening_period,closing_period)

            for j in xrange(opening_period,closing_period+1):
                nBlocksToExtract = self.units * schedule[i,j]
                blocks = dp.extraction(prevExtractions,nBlocksToExtract)
                extractedBlocks = len(blocks)
                prevExtractions += extractedBlocks

                #print "ton_blocks",i,j,self.units,schedule[i,j],extractedBlocks,np.sum(self.tonnage[blocks])


                if (extractedBlocks > 0):
                    '''tonnage'''
                    ton_blocks = np.sum(self.tonnage[blocks]) / 1000.0  #kilo tonns
                    #logger.debug("%d ton_blocks[%d-%d]=%f",nBlocksToExtract,i,j,ton_blocks)
                    
                    production_period[j] += ton_blocks


        return production_period


    def calculate_feed_tonnage_deviation(self,schedule):
        production_period = np.zeros(self.nperiods)

        for i in xrange(self.ndp):
            prevExtractions = 0

            dp = self.drawpoints[i]

            #logger.debug("processing DP: %d",i)

            for j in xrange(0,self.nperiods):
                nBlocksToExtract = self.units * schedule[i,j]
                blocks = dp.extraction(prevExtractions,nBlocksToExtract)
                extractedBlocks = len(blocks)
                prevExtractions += extractedBlocks

                if (extractedBlocks > 0):
                    '''tonnage'''
                    ton_blocks = np.sum(self.tonnage[blocks]) #/ 1000.0  #kilo tonns
                    #logger.debug("%d ton_blocks[%d-%d]=%f",nBlocksToExtract,i,j,ton_blocks)
                    
                    production_period[j] += ton_blocks

        #now compute deviation from min and max
        indices_less = np.where(production_period < self.minimum_feed_production)[0]
        indices_more = np.where(production_period > self.maximum_feed_production)[0]
        
        return production_period


    def compute_sim_objectives(self,schedule,nsim):
        #limit the computation for a one simulation at once
        nsr = self.nsr[:,nsim,:]
        return self.compute_nosim_objectives(schedule,nsr)

    def compute_nosim_objectives(self,schedule,nsr):
        assert(schedule.shape == (self.ndp,self.nperiods))

        npv_average = 0
        npv_std = 0
        grade_average = 0

        productionPeriod = np.zeros(self.nperiods)
        productionPeriod.fill(0.0);


        npv_sim = 0.0

        #blockExtracted = set()

        for i in xrange(self.ndp):
            prevExtractions = 0

            dp = self.drawpoints[i]

            #//cout << "processing DP: " << i<< endl;

            totalTonnage = 0.0

            for j in xrange(self.nperiods):

                nBlocksToExtract = self.units * schedule[i,j]
                 
                blocks = dp.extraction(prevExtractions,nBlocksToExtract)
                extractedBlocks = len(blocks)
                prevExtractions += extractedBlocks

                #// cout << "processing period: " << j<< endl;

                #print i,j,extractedBlocks
                
                if (extractedBlocks > 0):

                    nsr_blocks = nsr[:,blocks]
                    ton_blocks = self.tonnage[blocks]
                    npv_sum = np.sum((nsr_blocks * ton_blocks) * (self.discount[j] ))
                    
                    #print nsr_blocks.shape
                    #print ton_blokcs.shape
                    #print npv_sum.shape

                    npv_sim += npv_sum
                    
                    
                    productionPeriod[j] += np.sum(ton_blocks)


        dev = np.abs(productionPeriod - self.targetProduction)

        return npv_sim,np.sum(dev)

    def compute_nosim_openclose_objectives(self,schedule,opening,closing,nsr=None):
        '''
        this optimization method include opening and closing drawpoints
        '''
        assert(schedule.shape == (self.ndp,self.nperiods))
        assert(schedule.shape[0] == opening.shape[0])
        assert(schedule.shape[0] == closing.shape[0])
        
        if nsr is None:
            logger.debug("self.nsr_average.shape=%s",self.nsr_average.shape)
            nsr = self.nsr_average[:]

        npv_average = 0
        npv_std = 0
        grade_average = 0

        productionPeriod = np.zeros(self.nperiods)
        productionPeriod.fill(0.0);


        npv_sim = 0.0

        #blockExtracted = set()

        for i in xrange(self.ndp):
            prevExtractions = 0

            dp = self.drawpoints[i]

            #//cout << "processing DP: " << i<< endl;

            totalTonnage = 0.0

            opening_period = opening[i]
            closing_period = closing[i]

            #print i,opening_period,closing_period

            for j in xrange(opening_period,closing_period+1):

                nBlocksToExtract = self.units * schedule[i,j]
                
                blocks = dp.extraction(prevExtractions,nBlocksToExtract)
                extractedBlocks = len(blocks)
                prevExtractions += extractedBlocks

                #// cout << "processing period: " << j<< endl;

                #print i,j,extractedBlocks
                
                if (extractedBlocks > 0):

                    nsr_blocks = nsr[blocks]
                    ton_blocks = self.tonnage[blocks]
                    npv_sum = np.sum((nsr_blocks * ton_blocks) * (self.discount[j]))
                    
                    #print nsr_blocks.shape
                    #print ton_blokcs.shape
                    #print npv_sum.shape

                    npv_sim += npv_sum
                    
                    
                    productionPeriod[j] += np.sum(ton_blocks)


        dev = np.abs(productionPeriod - self.targetProduction)

        return (npv_sim,) #np.mean(dev)
