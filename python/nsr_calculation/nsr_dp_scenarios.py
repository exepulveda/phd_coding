import sys

sys.path += ["../base"]

import argparse
import random
import numpy as np
import geometry
import os.path
import pickle
import h5py
import matplotlib.pyplot as plt

def evaluate(individual):
    schedule = np.array(individual.reshape(bcp.ndp,bcp.nperiods),dtype=np.int)
    
    ret = bcp.calculate_maximum_npv(schedule,nsr_data)
    
    individual = schedule.flatten()
    
    return ret
    
parser = argparse.ArgumentParser()
parser.add_argument('--ipath', required=False,type=str,default="../results")    
parser.add_argument('--param', required=False,type=str,default="../data/pc1s1.json")    
    
if __name__ == "__main__":
    args = parser.parse_args()
    
    
    ds_path = args.ipath
    json_parameters = args.param

    if os.path.exists(os.path.join(ds_path,"grid3d.dump")):
        with open(os.path.join(ds_path,"grid3d.dump"),"r") as fin: 
            blockmodel3D = pickle.load(fin)


    #load dp
    dp_filename = os.path.join(ds_path,"drawpoints.mat")
    h5_dp = h5py.File(dp_filename, "r")
    dp_blocks_data = h5_dp["blockindices"][:,:]
    ndp,max_blocks = dp_blocks_data.shape
    h5_dp.close()

    #load tonnage
    density_filename = os.path.join(ds_path,"density.mat")
    h5_density = h5py.File(density_filename, "r")
    ds_tonnage = h5_density["density"][:] * blockmodel3D.volume()
    h5_density.close()
    
    
    
    nsr_filename = os.path.join(ds_path,"nsr_average_inputs.mat")
    h5_nsr = h5py.File(nsr_filename, "r")
    #load simulation
    nsr_data = h5_nsr["/average_input/nsr"][:]
    n,nsim = nsr_data.shape
    h5_nsr.close()

    dp_blocks = []
    for i in xrange(ndp):
        indices = np.where(dp_blocks_data[i,:] >= 0)[0]
        dp_blocks += [dp_blocks_data[i,indices]]

    nsr = np.zeros(nsim)
    tonnage = np.zeros(nsim)
    for i in xrange(nsim):
        n = 0
        for j in xrange(ndp):
            indices = dp_blocks[j]
            
            nsr[i] += np.sum(nsr_data[indices,i] * ds_tonnage[indices])
            tonnage[i] += np.sum(ds_tonnage[indices])
            n += len(indices)

        print (i+1),",",tonnage[i],",",nsr[i]
    #plt.hist(nsr,10)
    #plt.show()
