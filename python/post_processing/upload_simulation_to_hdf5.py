import sys

sys.path += ["../../../Documents/projects/geostatpy"]


#import decorrelating
import gaussian
import h5py
import numpy as np
import pickle
import os.path
import csv

if __name__ == "__main__":
    
    #load data for large number could be bad, better option is read and write by rows
    path = "/home/esepulveda/Documents/projects/newcrest/optimisation/ppmt"
    ds_path = "/home/esepulveda/Documents/projects/newcrest/optimisation/results"    
    
    reader = csv.reader(open(os.path.join(path,"ppmt_b.out"),"r"),delimiter=" ",skipinitialspace=True)
    #skip some gslib headers
    reader.next() #header
    row = reader.next() #variable number
    nvars = int(row[0])
    for _ in xrange(nvars):
        reader.next() #header
        
    if os.path.exists(os.path.join(ds_path,"grid3d.dump")):
        with open(os.path.join(ds_path,"grid3d.dump"),"r") as fin: 
            blockmodel3D = pickle.load(fin)

    batchsize = len(blockmodel3D)
    
    grade_simulations_filename = os.path.join(ds_path,"grade_simulations.mat")
    h5_grade_simulations = h5py.File(grade_simulations_filename, "r+")
    ds_au = h5_grade_simulations["/simulations/au"]
    ds_cu = h5_grade_simulations["/simulations/cu"]
    ds_s = h5_grade_simulations["/simulations/s"]
    ds_cucn = h5_grade_simulations["/simulations/cucn"]
    ds_fe = h5_grade_simulations["/simulations/fe"]
    ds_mo = h5_grade_simulations["/simulations/mo"]
    #au_ppm,cu_pct,s_pct,cucn_ppm,fe_pct,mo_ppm

    batch_data = np.empty((batchsize,6))
    nsim = 0
    k = 0
    for i,row in enumerate(reader):
        if i % batchsize == 0 and i > 0:
            print "current simulation",nsim
            ds_au[:,nsim] = batch_data[:,0]
            ds_cu[:,nsim] = batch_data[:,1]
            ds_s[:,nsim] = batch_data[:,2]
            ds_cucn[:,nsim] = batch_data[:,3]
            ds_fe[:,nsim] = batch_data[:,4]
            ds_mo[:,nsim] = batch_data[:,5]

            nsim += 1
            k = 0
            
        batch_data[k,:] = np.array([float(x) for x in row[:nvars]])
        k += 1

    #last sim
    print "current simulation",nsim
    ds_au[:,nsim] = batch_data[:,0]
    ds_cu[:,nsim] = batch_data[:,1]
    ds_s[:,nsim] = batch_data[:,2]
    ds_cucn[:,nsim] = batch_data[:,3]
    ds_fe[:,nsim] = batch_data[:,4]
    ds_mo[:,nsim] = batch_data[:,5]
    
