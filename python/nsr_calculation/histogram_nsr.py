import numpy as np
import sys
import os.path
import matplotlib.pyplot as plt
import h5py

ds_path = "/home/esepulveda/Documents/projects/newcrest/optimisation/nsr_results"

nsr_filename = os.path.join(ds_path,"nsr_A.mat")
h5nsr = h5py.File(nsr_filename, "r")

dset_nsr = h5nsr["nsr"]

if len(dset_nsr.shape) == 1:
    data = dset_nsr[:]
    print np.min(data)
    print np.max(data)
    print np.mean(data)    
    data = np.clip(data,-100,100)
    
    
    plt.hist(data,11,normed=True)
    plt.show()

quit()    

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
