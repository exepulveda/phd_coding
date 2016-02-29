import numpy as np
import sys
import os.path
import matplotlib.pyplot as plt
import h5py

ds_path = "/home/esepulveda/Documents/projects/newcrest/optimisation/results"

nsr_filename = os.path.join(ds_path,"nsr_average_inputs.mat")
h5nsr = h5py.File(nsr_filename, "r")

dset_nsr = h5nsr["average_input/nsr"]
dset_nsr_average = h5nsr["/average_input/nsr_average"]

n,m = dset_nsr.shape

for i in xrange(m):
    data = dset_nsr[:,i]
    #plt.hist(data,9,normed=True)
    #plt.show()
    #print (i+1),np.min(data),np.mean(data),np.max(data)
    

for i in xrange(n):
    data = dset_nsr[i,:]
    #plt.hist(data,9,normed=True)
    #plt.show()
    print (i+1),np.min(data),np.mean(data),np.max(data)
