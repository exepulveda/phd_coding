import sys

sys.path += ["/home/esepulveda/Documents/projects/geostatpy"]
sys.path += ["../base"]

import argparse
import random
import numpy as np
import multiprocessing
import geometry
import os.path
import pickle
import h5py
import bcproblem

def cxTwoPointCopy(ind1, ind2):
    """Execute a two points crossover with copy on the input individuals. The
    copy is required because the slicing in numpy returns a view of the data,
    which leads to a self overwritting in the swap operation. It prevents
    ::
    
        >>> import numpy
        >>> a = numpy.array((1,2,3,4))
        >>> b = numpy.array((5.6.7.8))
        >>> a[1:3], b[1:3] = b[1:3], a[1:3]
        >>> print(a)
        [1 6 7 4]
        >>> print(b)
        [5 6 7 8]
    """
    size = len(ind1)
    cxpoint1 = random.randint(1, size)
    cxpoint2 = random.randint(1, size - 1)
    if cxpoint2 >= cxpoint1:
        cxpoint2 += 1
    else: # Swap the two cx points
        cxpoint1, cxpoint2 = cxpoint2, cxpoint1

    ind1[cxpoint1:cxpoint2], ind2[cxpoint1:cxpoint2] \
        = ind2[cxpoint1:cxpoint2].copy(), ind1[cxpoint1:cxpoint2].copy()
        
    return ind1, ind2


def initScheduler():
    return np.random.randint(0,200)

def evaluate(individual):
    schedule = np.array(individual.reshape(bcp.ndp,bcp.nperiods),dtype=np.int)
    
    ret = bcp.calculate_maximum_npv(schedule,nsr_data)
    
    individual = schedule.flatten()
    
    return ret
    
parser = argparse.ArgumentParser()
parser.add_argument('--periods', required=False,type=int,default=12)    
parser.add_argument('--pop', required=False,type=int,default=100)    
parser.add_argument('--gen', required=False,type=int,default=100)    
parser.add_argument('--cxpb', required=False,type=float,default=0.8)    
parser.add_argument('--mutpb', required=False,type=float,default=0.2)    
parser.add_argument('--ipath', required=False,type=str,default="../results")    
parser.add_argument('--param', required=False,type=str,default="../data/pc1s1.json")    
parser.add_argument('--seed', required=False,type=int,default=1634120)    
    
if __name__ == "__main__":
    args = parser.parse_args()
    
    random.seed(args.seed)
    np.random.seed(args.seed)
    
    output_path = "./output-optimisation"

    ds_path = args.ipath

    json_parameters = args.param

    from bcproblem import BlockCavingProblem
        
    bcp = BlockCavingProblem()
    bcp.load(json_parameters,load_nsr=False)

    if os.path.exists(os.path.join(ds_path,"grid3d.dump")):
        with open(os.path.join(ds_path,"grid3d.dump"),"r") as fin: 
            blockmodel3D = pickle.load(fin)


    #load dp
    dp_filename = os.path.join(ds_path,"drawpoints.mat")
    h5_dp = h5py.File(dp_filename, "r")
    dp_blocks_data = h5_dp["blockindices"][:,:]
    ndp,max_blocks = dp_blocks_data.shape
    h5_dp.close()

    nsr_filename = os.path.join(ds_path,"best_nsr_average_inputs.mat")
    h5_nsr = h5py.File(nsr_filename, "r")
    #load simulation
    nsr_data = h5_nsr["/average_input/nsr"][:]
    h5_nsr.close()


    
    dp_blocks = []
    for i in xrange(ndp):
        indices = np.where(dp_blocks_data[i,:] >= 0)[0]
        dp_blocks += [dp_blocks_data[i,indices]]

    nperiods = args.periods

    bcp.setup_drawpoints(dp_blocks)
    bcp.setup_periods(nperiods,0.1)

    #schedule = np.random.randint(0,50,size=(ndp,nperiods))
    
    #ret = bcp.compute_objectives(schedule)
    #print ret

    from deap import tools
    from deap import base, creator
    from deap import algorithms

    creator.create("FitnessNSR", base.Fitness, weights=(1.0,))
    creator.create("Individual", np.ndarray, fitness=creator.FitnessNSR)

    toolbox = base.Toolbox()

    #pool = multiprocessing.Pool()
    #toolbox.register("map", pool.map)


    #toolbox.register("map", futures.map)

    toolbox.register("mate", cxTwoPointCopy)
    #toolbox.register("mutate", tools.mutGaussian, mu=0, sigma=10, indpb=0.1)
    toolbox.register("mutate", tools.mutPolynomialBounded,eta=1,low=0, up=200, indpb=0.05)
    toolbox.register("select", tools.selTournament, tournsize=5)
    #toolbox.register("select", tools.selRoulette)
    toolbox.register("evaluate", evaluate)
    

    IND_SIZE = bcp.ndp * bcp.nperiods

    toolbox.register("attribute", initScheduler)
    toolbox.register("individual", tools.initRepeat, creator.Individual,
                     toolbox.attribute, n=IND_SIZE)
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)

    NGEN = args.gen
    MU = args.pop
    CXPB = args.cxpb
    MUTPB = args.mutpb
    
    pop = toolbox.population(n=MU)
    hof = tools.HallOfFame(10,similar=np.array_equal)
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", np.mean, axis=0)
    stats.register("std", np.std, axis=0)
    stats.register("min", np.min, axis=0)
    stats.register("max", np.max, axis=0)
    
    #stats = None
    print "Start evolution!!"
    algorithms.eaSimple(pop, toolbox, cxpb=CXPB, mutpb=MUTPB, ngen=NGEN, stats=stats,
                              halloffame=hof,verbose=True)

    #generate output
    try:
        os.makedirs(output_path)
    except:
        pass
    
    for i,ind in enumerate(hof):
        solution_path = os.path.join(output_path,"solution-{0}.csv".format(i+1))
        np.savetxt(solution_path,ind.reshape(bcp.ndp,bcp.nperiods),fmt="%d",delimiter=",")
        
        print i+1,ind.fitness.values[0]
