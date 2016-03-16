import sys

sys.path += ["/home/esepulveda/Documents/projects/geostatpy"]
sys.path += ["../base"]

import argparse
import random
import numpy as np
import geometry
import os.path
import pickle
import h5py
import bcproblem

parser = argparse.ArgumentParser()
parser.add_argument('--solution_home', required=True,type=str)    
parser.add_argument('--solution_name', required=False,type=str,default="output-geomet-sim-{0}")    
parser.add_argument('--ipath', required=True,type=str)    
parser.add_argument('--param', required=False,type=str,default="pc1s1.json")    
parser.add_argument('--nsr_file', required=False,type=str,default="nsr_B.mat")    
parser.add_argument('--nsr_ds', required=False,type=str,default="nsr")    
    
if __name__ == "__main__":
    args = parser.parse_args()
    
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

    nsr_filename = os.path.join(ds_path,args.nsr_file)
    h5_nsr = h5py.File(nsr_filename, "r")
    #load simulation
    #h5_nsr.close()

    dp_blocks = []
    for i in xrange(ndp):
        indices = np.where(dp_blocks_data[i,:] >= 0)[0]
        dp_blocks += [dp_blocks_data[i,indices]]

    nsim = 50
    #load solutions
    for i in xrange(nsim):
        solution_file = os.path.join(args.solution_home,args.solution_name.format(i),"solution-1.csv")
        #print "processing",solution_file,"sim",i
        schedule = np.loadtxt(solution_file,delimiter=',')
        assert (ndp == schedule.shape[0])
        
        _,nperiods = schedule.shape

        bcp.setup_drawpoints(dp_blocks)
        bcp.setup_periods(nperiods,0.1)

        nsr_data = h5_nsr[args.nsr_ds][:,i]


        npv,tonnage = bcp.calculate_npv_tonnage(schedule,nsr_data)
        
        print (i+1),",",np.sum(npv)            
        
