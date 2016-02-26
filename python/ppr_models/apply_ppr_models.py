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
    ("au_recleaner_concentrate","au/concentrated"),
    ("cu_recleaner_concentrate","cu/concentrated"),
    ("f_recleaner_concentrate","f/concentrated"),
    ("au_rougher_recovery","au/recovery"),
    ("cu_rougher_recovery","cu/recovery")  ]

clip_values = {
    "au_recleaner_concentrate":(3.0,),
    "cu_recleaner_concentrate":(10.0,),
    "f_recleaner_concentrate":(100.0,),
    "au_rougher_recovery":(30.0,),
    "cu_rougher_recovery":(40.0,)
    ]

ngeomet = 50

ppr_path = "/home/esepulveda/Documents/projects/newcrest/optimisation/ppr_models"

ppr_model = {}
for modelname,_ in modelnames:
    ppr_model[modelname] = []
    for i in xrange(ngeomet):
        fn = os.path.join(ppr_path,"ppr_model_{name}_{number}.dump")
        with open(fn.format(name=modelname,number=i),"r") as fout:
            ppr_model[modelname] += [pickle.load(fout)]

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
nsr_filename = os.path.join(ds_path,"nsr_average_inputs.mat")
#h5nsr = h5py.File(nsr_filename, "r+")
#for modelname,ds_name in modelnames:
#    h5nsr.create_dataset(ds_name, (n,ngeomet), dtype=np.float)
#
#h5nsr.close()

h5nsr = h5py.File(nsr_filename, "r+")
for modelname,ds_name in modelnames:
    print "processing",ds_name
    ds = h5nsr[ds_name]
    for i in xrange(ngeomet):
        print "processing model nro",(i+1)
        ret = ppr_model[modelname][i].predict(inputs)
        
        ds[:,i] = ret
                
quit()
#simulation dataset

dset_aurec = h5_simulations["/au/recovery"]
n,m = dset_aurec.shape
print n,m

for j in xrange(m):
    print "realisation",j,
    print np.min(dset_aurec[:,j]),
    print np.mean(dset_aurec[:,j]),
    print np.max(dset_aurec[:,j])
