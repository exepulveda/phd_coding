'''configure datasource for optimisation
'''
import sys

sys.path += ["/home/esepulveda/Documents/projects/geostatpy"]
sys.path += ["/home/esepulveda/Documents/projects/newcrest/scripts/python"]
sys.path += ["/home/esepulveda/Documents/projects/fortran"]

import numpy as np
import scipy.spatial
import os.path
import csv
import pickle
import geometry
import h5py

ds_path = "/home/esepulveda/Documents/projects/newcrest/optimisation/results"

reload_density = False

#define or load grid
blockmodel3D = None
if os.path.exists(os.path.join(ds_path,"grid3d.dump")):
    with open(os.path.join(ds_path,"grid3d.dump"),"r") as fin: 
        try:
            blockmodel3D = pickle.load(fin)
            print "block model definition loaded"
        except:
            blockmodel3D = None
            
modelnames = [
    ("au_recleaner_concentrate","au/concentrated"),
    ("cu_recleaner_concentrate","cu/concentrated"),
    ("f_recleaner_concentrate","f/concentrated"),
    ("au_rougher_recovery","au/recovery"),
    ("cu_rougher_recovery","cu/recovery")  ]

clip_values = {
    "au_recleaner_concentrate":(3.0,600.0),
    "cu_recleaner_concentrate":(10.0,60.0),
    "f_recleaner_concentrate":(100.0,3000),
    "au_rougher_recovery":(40.0,95.0),
    "cu_rougher_recovery":(50.0,95.0)
    }

nsr_filename = os.path.join(ds_path,"nsr_average_inputs.mat")
h5nsr = h5py.File(nsr_filename, "r+")

for modelname,ds_name in modelnames:
    print "processing",ds_name
    ds = h5nsr[ds_name]
    n,m = ds.shape
    dmin,dmax = clip_values[modelname]
    for i in xrange(m):
        data = ds[:,i]
        print "processing model nro",(i+1),np.min(data),np.mean(data),np.max(data)
        data = np.clip(data,dmin,dmax)
        ds[:,i] = data
        print "processing model nro",(i+1),np.min(data),np.mean(data),np.max(data)
