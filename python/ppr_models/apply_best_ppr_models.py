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
grade_simulations_filename = os.path.join(ds_path,"average_grade_simulations.mat")
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

ppr_path = "/home/esepulveda/Documents/projects/newcrest/optimisation/ppr_models"

ppr_model = {}
for modelname,_,i in modelnames:
    ppr_model[modelname] = []
    fn = os.path.join(ppr_path,"ppr_model_{name}_{number}.dump")
    with open(fn.format(name=modelname,number=i),"r") as fout:
        ppr_model[modelname] = pickle.load(fout)

ds_inputs = h5_grade_simulations["/simulations"]

n = len(ds_inputs["au"])

assert (len(blockmodel3D) == n)

inputs = np.empty((7,n))

au_inputs = ds_inputs["au"][:]
cu_inputs = ds_inputs["cu"][:]
s_inputs = ds_inputs["s"][:]
cucn_inputs = ds_inputs["cucn"][:]
fe_inputs = ds_inputs["fe"][:]
mo_inputs = ds_inputs["mo"][:]
f_inputs = ds_inputs["f"][:]

inputs.fill(-999)

inputs[0,:] = au_inputs
inputs[1,:] = cu_inputs
inputs[2,:] = mo_inputs
inputs[3,:] = fe_inputs
inputs[4,:] = s_inputs
inputs[5,:] = f_inputs
inputs[6,:] = cucn_inputs / 10000.0
    
    #    inputs = [au_ppm,cu_pct,mo_ppm,fe_pct,s_pct,f_ppm,cucn_pct]

#create and populate nsr models with average inputs
nsr_filename = os.path.join(ds_path,"best_nsr_average_inputs.mat")
h5nsr = h5py.File(nsr_filename, "r+")
for modelname,ds_name,i in modelnames:
    dmin,dmax = clip_values[modelname]
    print "processing",ds_name,(i+1),"clipping",dmin,dmax
    #ds = h5nsr.create_dataset(ds_name, (n,), dtype='f')
    #ds = h5nsr[ds_name]
    ret = ppr_model[modelname].predict(inputs)
    ret = np.clip(ret,dmin,dmax)
    ds[:] = ret
                
