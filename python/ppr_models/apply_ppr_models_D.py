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


clip_values = {
    "au_recleaner_concentrate":(3.0,600.0),
    "cu_recleaner_concentrate":(10.0,60.0),
    "f_recleaner_concentrate":(100.0,3000),
    "au_rougher_recovery":(40.0,95.0),
    "cu_rougher_recovery":(50.0,95.0)
    }

ds_inputs = h5_grade_simulations["/simulations"]
n,nsim = ds_inputs["au"].shape

#nsim = 
ngeomet = 50

print "simulations shape",n,nsim,ngeomet

ppr_path = "/home/esepulveda/Documents/projects/newcrest/optimisation/ppr_models"
ppr_model = {}
for modelname,_,_ in modelnames:
    ppr_model[modelname] = []
    for i in xrange(ngeomet):
        fn = os.path.join(ppr_path,"ppr_model_{name}_{number}.dump")
        with open(fn.format(name=modelname,number=i),"r") as fout:
            ppr_model[modelname] += [pickle.load(fout)]

assert (len(blockmodel3D) == n)

au_inputs = ds_inputs["au"][:,:]
cu_inputs = ds_inputs["cu"][:,:]
s_inputs = ds_inputs["s"][:,:]
cucn_inputs = ds_inputs["cucn"][:,:]
fe_inputs = ds_inputs["fe"][:,:]
mo_inputs = ds_inputs["mo"][:,:]
f_inputs = ds_inputs["f"][:,:]

inputs = np.empty((7,n))


#create and populate geomet models with 50 grades and 50 geomet simulations
nsr_filename = os.path.join(ds_path,"geomet_D.mat")
if True:
    h5nsr = h5py.File(nsr_filename, "a")
    for modelname,ds_name,_ in modelnames:
        h5nsr.create_dataset(ds_name, (n,nsim*ngeomet), dtype=np.float)

    h5nsr.close()


h5nsr = h5py.File(nsr_filename, "r+")

for k in xrange(nsim):
    inputs.fill(-999)

    inputs[0,:] = au_inputs[:,k]
    inputs[1,:] = cu_inputs[:,k]
    inputs[2,:] = mo_inputs[:,k]
    inputs[3,:] = fe_inputs[:,k]
    inputs[4,:] = s_inputs[:,k]
    inputs[5,:] = f_inputs[:,k]
    inputs[6,:] = cucn_inputs[:,k] / 10000.0
    
    print "processing simulation",(k+1),"..."

    for modelname,ds_name,_ in modelnames:
        print "processing",ds_name,ngeomet
        ds = h5nsr[ds_name]
        dmin,dmax = clip_values[modelname]
        for i in xrange(ngeomet):
            ret = ppr_model[modelname][i].predict(inputs)
            
            simulation = k*ngeomet + i
            
            ds[:,simulation] = np.clip(ret,dmin,dmax)
                
