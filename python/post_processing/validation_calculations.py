import sys

sys.path += ["/home/esepulveda/Documents/projects/geostatpy"]

import random
import numpy as np
import multiprocessing
import logging
import os.path
import pickle
import geometry
import h5py
import csv
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    #logging.basicConfig(level=logging.DEBUG)

    json_parameters = "/home/esepulveda/Documents/projects/newcrest/optimisation/data/pc1s1.json";

    path = "/home/esepulveda/Documents/projects/newcrest/optimisation/sgsim/results"
    ds_path = "/home/esepulveda/Documents/projects/newcrest/optimisation/results"   
    
    #filename = sys.argv[1]

    if os.path.exists(os.path.join(ds_path,"grid3d.dump")):
        with open(os.path.join(ds_path,"grid3d.dump"),"r") as fin: 
            blockmodel3D = pickle.load(fin)

    dp_filename = os.path.join(ds_path,"drawpoints.mat")
    h5_dp = h5py.File(dp_filename, "r")

    dp_filename = os.path.join(ds_path,"grade_simulations.mat")
    h5_grades = h5py.File(dp_filename, "r")
    
    for var in ["au","cu","f","fe","mo","s","cucn"]:
        ds_name = "simulations/" + var
        ds = h5_grades[ds_name][:,:]
        n,nsim = ds.shape

        #for sim in xrange(nsim):
        #    print var,sim,np.min(ds[:,sim]),np.mean(ds[:,sim]),np.max(ds[:,sim]),np.std(ds[:,sim])

        #print var,sim,np.min(ds[:,:]),np.mean(ds[:,:]),np.max(ds[:,:]),np.std(ds[:,:])

        #plt.hist(ds,bins=31,normed=True)
        #plt.title(var)
        #plt.show()

        #ds_std = np.std(ds,axis=1)
        #print var,np.mean(ds_std)
    
        index = 45000
        data = ds[index,:]
        print var,"mean",np.mean(data),"std",np.std(data)
        plt.hist(data,bins=9,normed=True)
        plt.title(var)
        plt.show()

        #quit()
