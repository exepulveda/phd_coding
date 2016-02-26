'''create average simulations values
'''
import sys
import numpy as np
import os.path
import csv
import h5py

ds_path = "/home/esepulveda/Documents/projects/newcrest/optimisation/results"

#open grade_simulations
simulations_filename = os.path.join(ds_path,"grade_simulations.mat")
h5_simulations = h5py.File(simulations_filename, "r")

#dataosurces
ds_names = ["au","cu","s","cucn","fe","mo","f"]

#create grade_simulation dataset
average_grade_simulations_filename = os.path.join(ds_path,"average_grade_simulations.mat")
h5_average_grade_simulations = h5py.File(average_grade_simulations_filename, "w")

for ds_name in ds_names:
    data = h5_simulations["/simulations/" + ds_name][:,:]
    n,m = data.shape
    print "processing",ds_name,n,m
    ds = h5_average_grade_simulations.create_dataset("/simulations/" + ds_name, (n,), dtype=np.float)
    ds[:] = np.mean(data,axis=1)
