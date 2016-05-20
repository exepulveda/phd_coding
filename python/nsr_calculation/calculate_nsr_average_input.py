import numpy as np
import sys
import h5py
import json
import os.path
import pickle

sys.path += ["/home/esepulveda/Documents/projects/geostatpy"]
sys.path += ["/home/esepulveda/Documents/projects/newcrest/scripts/python"]
sys.path += ["/home/esepulveda/Documents/projects/fortran"]

import csv
import geometry
import h5py

import nsr

ds_path = "/home/esepulveda/Documents/projects/newcrest/optimisation/results"

#define or load grid
blockmodel3D = None
if os.path.exists(os.path.join(ds_path,"grid3d.dump")):
    with open(os.path.join(ds_path,"grid3d.dump"),"r") as fin: 
        try:
            blockmodel3D = pickle.load(fin)
            print "block model definition loaded"
        except Exception as e:
            print "block model definition loaded",e
            blockmodel3D = None

nx,ny,nz = blockmodel3D.nodes

ds_path = "/home/esepulveda/Documents/projects/newcrest/optimisation/nsr_results"

nblocks = nx*ny*nz

#open nsr sim
nsr_filename = os.path.join(ds_path,"nsr_D.mat")
h5nsr = h5py.File(nsr_filename, "r")

#create nsr dataset
nsr_average_filename = os.path.join(ds_path,"nsr_D_average.mat")
h5nsr_average = h5py.File(nsr_average_filename, "w")

if True:
    dset_nsr_average = h5nsr_average.create_dataset("nsr", (nblocks,), dtype='f')

dset_nsr = h5nsr["nsr"]
n,nsim = dset_nsr.shape
dset_nsr_average = h5nsr_average["nsr"]

dset_nsr_average[:] = np.mean(dset_nsr[:,:],axis=1)

#split
for i in range(nsim):
    nsr_filename = os.path.join(ds_path,"nsr_D_{0}.mat".format(i))
    h5nsr_sim = h5py.File(nsr_filename, "w")
    dset_nsr_sim = h5nsr_sim.create_dataset("nsr", (nblocks,), dtype='f')
    dset_nsr_sim[:] = dset_nsr[:,i]
    h5nsr_sim.close()
    
    
