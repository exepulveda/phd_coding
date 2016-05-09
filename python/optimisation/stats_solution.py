import sys

sys.path += ["/home/esepulveda/Documents/projects/geostatpy"]
sys.path += ["/home/esepulveda/phd_coding/python/base"]

import argparse
import random
import numpy as np
import geometry
import os.path
import pickle
import h5py
import bcproblem

parser = argparse.ArgumentParser()
parser.add_argument('--solution_file', required=True,type=str)    
parser.add_argument('--ipath', required=True,type=str)    
parser.add_argument('--param', required=False,type=str,default="pc1s1.json")    
parser.add_argument('--nsr_file', required=False,type=str,default="best_nsr_average_inputs.mat")    
parser.add_argument('--nsr_ds', required=False,type=str,default="/average_input/nsr")
parser.add_argument('--sim', required=False,type=int,default=-1)
    
if __name__ == "__main__":
    example = '''
        python stats_solution.py --solution_file=.mat --param=pc1s1_A.json --nsr_file= --nsr_ds=nsr 
    '''
    
    
    args = parser.parse_args()
    
    output_path = "./output-optimisation"

    ds_path = args.ipath
    json_parameters = args.param

    from bcproblem import BlockCavingProblem
        
    bcp = BlockCavingProblem()
    bcp.load(json_parameters,load_nsr=False)

    if os.path.exists(os.path.join(ds_path,"grid3d.dump")):
        with open(os.path.join(ds_path,"grid3d.dump"),"r") as fin: 
            blockmodel3D = pickle.load(fin)

    nsr_filename = os.path.join(ds_path,args.nsr_file)
    h5_nsr = h5py.File(nsr_filename, "r")
    #load simulation
    if args.sim >= 0:
        nsr_data = h5_nsr[args.nsr_ds][:,args.sim]
    else:
        nsr_data = h5_nsr[args.nsr_ds][:]
    h5_nsr.close()

    #load solution
    schedule = np.loadtxt(args.solution_file,delimiter=',')
    
    ndp,nperiods = schedule.shape

    npv,tonnage = bcp.calculate_npv_tonnage(schedule,nsr_data)
    
    print "NSR:",np.sum(npv)
    print "Tonnage:",np.sum(tonnage)

    print "Period","Accumulated NSR","Tonnage"
    
    npv = np.cumsum(npv)
    for i in xrange(nperiods):
        print (i+1),npv[i],tonnage[i]
        
    
