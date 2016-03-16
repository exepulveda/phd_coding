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
            
#average grades
grade_simulations_filename = os.path.join(ds_path,"grade_simulations.mat")
h5_grade_simulations = h5py.File(grade_simulations_filename, "r")

modelnames = [
    ("au_recleaner_concentrate","au/concentrated",1),
    ("cu_recleaner_concentrate","cu/concentrated",35),
    ("f_recleaner_concentrate","f/concentrated",2),
    ("au_rougher_recovery","au/recovery",7),
    ("cu_rougher_recovery","cu/recovery",23)  ]

modelnames = [
    ("au_recleaner_concentrate","au/concentrated",6),  #1,35,2,7,23
    ("cu_recleaner_concentrate","cu/concentrated",3),
    ("f_recleaner_concentrate","f/concentrated",1),
    ("au_rougher_recovery","au/recovery",17),
    ("cu_rougher_recovery","cu/recovery",23)  ]

clip_values = {
    "au_recleaner_concentrate":(3.0,600.0),
    "cu_recleaner_concentrate":(10.0,60.0),
    "f_recleaner_concentrate":(100.0,3000),
    "au_rougher_recovery":(40.0,95.0),
    "cu_rougher_recovery":(50.0,95.0)
    }

ds_inputs = h5_grade_simulations["/simulations"]
n,nsim = ds_inputs["au"].shape

print "simulations shape",n,nsim

ppr_path = "/home/esepulveda/Documents/projects/newcrest/optimisation/ppr_models_5"
ppr_model = {}
for modelname,_, ngeomet in modelnames:
    fn = os.path.join(ppr_path,"ppr_model_{name}_{number}.dump")
    with open(fn.format(name=modelname,number=ngeomet),"r") as fout:
        ppr_model[modelname] = pickle.load(fout)

assert (len(blockmodel3D) == n)

au_inputs = ds_inputs["au"][:,:]
cu_inputs = ds_inputs["cu"][:,:]
s_inputs = ds_inputs["s"][:,:]
cucn_inputs = ds_inputs["cucn"][:,:]
fe_inputs = ds_inputs["fe"][:,:]
mo_inputs = ds_inputs["mo"][:,:]
f_inputs = ds_inputs["f"][:,:]

#inputs = np.empty((7,n))
inputs = np.empty((6,n))


#create and populate nsr models with average inputs
nsr_filename = os.path.join(ds_path,"geomet_B.mat")
if False:
    h5nsr = h5py.File(nsr_filename, "a")
    for modelname,ds_name,_ in modelnames:
        h5nsr.create_dataset(ds_name, (n,nsim), dtype=np.float)

    h5nsr.close()

h5nsr = h5py.File(nsr_filename, "r+")

ret_data = np.empty((n,nsim))

for modelname,ds_name,ngeomet in modelnames:
    print "processing",ds_name,ngeomet
    for k in xrange(nsim):
        inputs[0,:] = au_inputs[:,k]
        inputs[1,:] = cu_inputs[:,k]
        inputs[2,:] = mo_inputs[:,k] / 10000.0
        inputs[3,:] = fe_inputs[:,k]
        inputs[4,:] = s_inputs[:,k]
        #inputs[5,:] = f_inputs[:,k]
        #inputs[6,:] = cucn_inputs[:,k] / 10000.0
        inputs[5,:] = cucn_inputs[:,k] / 10000.0

        dmin,dmax = clip_values[modelname]
        ret = ppr_model[modelname].predict(inputs)
        ret_data[:,k] = np.clip(ret,dmin,dmax)

    ds = h5nsr[ds_name]
    ds[:,:] = ret_data
