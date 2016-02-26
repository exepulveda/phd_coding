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
            
#grade
grade_simulations_filename = os.path.join(ds_path,"grade_simulations.mat")
h5_grade_simulations = h5py.File(grade_simulations_filename, "r")

model_names = ["au_recleaner_concentrate","cu_recleaner_concentrate","f_recleaner_concentrate","au_rougher_recovery","cu_rougher_recovery"]

ngeomet = 10

ppr_model = {}
for modelname in model_names:
    ppr_model[modelname] = []
    for i in xrange(ngeomet):
        with open("../models/ppr_model_{name}_{number}.dump".format(name=modelname,number=i),"r") as fout:
            ppr_model[modelname] += [pickle.load(fout)]

ds_inputs = h5_grade_simulations["/simulations"]

nsim = 50

inputs = np.empty((7,nsim))

nlog = 1000
for k in xrange(100):
    au_inputs = ds_inputs["au"][k,:]
    cu_inputs = ds_inputs["cu"][k,:]
    s_inputs = ds_inputs["s"][k,:]
    cucn_inputs = ds_inputs["cncu"][k,:]
    fe_inputs = ds_inputs["fe"][k,:]
    mo_inputs = ds_inputs["mo"][k,:]
    f_inputs = ds_inputs["f"][k,:]

    modelname = "au_recleaner_concentrate"

    
    inputs.fill(-999)

    inputs[0,:] = au_inputs
    inputs[1,:] = cu_inputs
    inputs[2,:] = mo_inputs
    inputs[3,:] = fe_inputs
    inputs[4,:] = s_inputs
    inputs[5,:] = f_inputs
    inputs[6,:] = cucn_inputs / 10000.0
    
    #    inputs = [au_ppm,cu_pct,mo_ppm,fe_pct,s_pct,f_ppm,cucn_pct]
    
    for i in xrange(ngeomet):
        inputs_i = inputs[:,i]
        ret = ppr_model[modelname][0].predict(inputs_i.reshape((7,1)))
        if ret < 0:
            print inputs,
        print modelname,i,ret
    
    np.savetxt('inputs.csv',inputs.T,fmt="%.6f",delimiter=",")
    
                
quit()
#simulation dataset
simulations_filename = os.path.join(ds_path,"simulations.mat")
h5_simulations = h5py.File(simulations_filename, "r")

dset_aurec = h5_simulations["/au/recovery"]
n,m = dset_aurec.shape
print n,m

for j in xrange(m):
    print "realisation",j,
    print np.min(dset_aurec[:,j]),
    print np.mean(dset_aurec[:,j]),
    print np.max(dset_aurec[:,j])
