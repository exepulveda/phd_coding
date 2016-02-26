import numpy as np
import sys
import h5py
import json
import os.path
import pickle

sys.path += ["/home/esepulveda/Documents/projects/geostatpy"]
sys.path += ["/home/esepulveda/Documents/projects/newcrest/scripts/python"]
sys.path += ["/home/esepulveda/Documents/projects/fortran"]

import csv
import geometry
import h5py

import nsr

ds_path = "/home/esepulveda/Documents/projects/newcrest/optimisation/results"

reload_density = False

#define or load grid
blockmodel3D = None
if os.path.exists(os.path.join(ds_path,"grid3d.dump")):
    with open(os.path.join(ds_path,"grid3d.dump"),"r") as fin: 
        try:
            blockmodel3D = pickle.load(fin)
            print "block model definition loaded"
        except Exception as e:
            print "block model definition loaded",e
            blockmodel3D = None

nx,ny,nz = blockmodel3D.nodes

n_simulations = 50
n_geomet = 10
nblocks = nx*ny*nz

#create nsr dataset
nsr_filename = os.path.join(ds_path,"nsr.mat")
h5_nsr = h5py.File(nsr_filename, "a")

if False:
    h5_nsr.create_dataset("/nsr", (nblocks,n_simulations*n_geomet), dtype=np.float)
    h5_nsr.create_dataset("/nsr_average", (nblocks,), dtype=np.float)

dset_nsr = h5_nsr["/nsr"]
dset_nsr_average = h5_nsr["/nsr_average"]

home = "/home/esepulveda/Documents/projects/newcrest/optimisation"

json_parameters = os.path.join(home,"data","pc1s1.json");
root = json.load(open(json_parameters))

h5grades = h5py.File(os.path.join(home,"results","grade_simulations.mat"),"a")
h5simulations = h5py.File(os.path.join(home,"results","simulations.mat"),"a")
h5nsr = h5py.File(os.path.join(home,"results","nsr.mat"),"a")

#cu datasets
dset_grades_cu = h5grades["simulations/cu"]
dset_concentrates_cu = h5simulations["cu/concentrated"]
dset_recovery_cu = h5simulations["cu/recovery"]

#au datasets
dset_grades_au = h5grades["simulations/au"]
dset_concentrates_au = h5simulations["au/concentrated"]
dset_recovery_au = h5simulations["au/recovery"]

#fluorine
dset_concentrates_f = h5simulations["f/concentrated"]

#batch size
def generate_batch(n,bacth_size):
    start = 0
    ret = []
    while start < n:
        ret += [(start, min(start+bacth_size,n))]
        start += bacth_size
        
    return ret


metal_node = root["metals"]["cu"]
price_cu = metal_node["price"]
price_to_ton_cu = metal_node["unit2ton"]
deduct_cu = metal_node["deduction"]
deduct_proportion_cu = metal_node["payableRate"]
refining_cost_cu = metal_node["refiningCharge"] #US$/lb refiningCharge

metal_node = root["metals"]["au"]
price_au = metal_node["price"]
price_to_ton_au = metal_node["unit2ton"]
deduct_au = metal_node["deduction"]
deduct_proportion_au = metal_node["payableRate"]
refining_cost_au = metal_node["refiningCharge"] #US$/oz

smelter_charges = root["smelterCharge"]
freight_charges = root["freightCharge"]
penalties = root["penalties"]
verbose = False


bacth_size = 10000

for s,e in generate_batch(nblocks,bacth_size):
    print "processing",s,e

    grade_cu = np.empty((e-s,n_simulations*n_geomet))
    grade_au = np.empty((e-s,n_simulations*n_geomet))
    
    grades_cu_tmp = dset_grades_cu[s:e,:]
    grades_au_tmp = dset_grades_au[s:e,:]
    concentrate_cu = dset_concentrates_cu[s:e,:]
    concentrate_au = dset_concentrates_au[s:e,:]
    recovery_cu = dset_recovery_cu[s:e,:]
    recovery_au = dset_recovery_au[s:e,:]
    
    #broadcast grades
    for i in xrange(n_geomet):
        grade_cu[:,i*n_simulations:(i+1)*n_simulations] = grades_cu_tmp[:,:]
        grade_au[:,i*n_simulations:(i+1)*n_simulations] = grades_au_tmp[:,:]

    del grades_cu_tmp
    del grades_au_tmp
    
        #marketing_charge = root["marketing"]

    nsr_ore_cu,nsr_concentrate_cu,nsr_ore_au,nsr_concentrate_au = nsr.calculateNSRLane(
        #cu
        grade_cu,
        recovery_cu,
        concentrate_cu,
        price_cu,
        price_to_ton_cu,
        #au
        grade_au,
        recovery_au,
        concentrate_au,
        price_au,
        price_to_ton_au,
        #cost cu
        deduct_cu,
        deduct_proportion_cu,
        refining_cost_cu,
        #cost au
        deduct_au,
        deduct_proportion_au,
        refining_cost_au,
        #common costs
        smelter_charges,
        freight_charges,
        verbose = True)        
    
    dset_nsr[s:e,:] = nsr_ore_cu + nsr_ore_au
    dset_nsr_average[s:e] = np.mean(nsr_ore_cu + nsr_ore_au,axis=1)
    
    
quit()

with h5py.File("../data/bm_{0}_{1}_{2}_dataset.hdf5".format(nblocks,n_simulations,n_geomet), "a") as f:
    #Cu
    #grades
    dset_grades_cu = f["simulation_cu"][:,:]
    #concentrates
    dset_concentrates_cu = f["geomet_cu_con"][:,:,:]
    #recovery
    dset_recovery_cu = f["geomet_cu_rec"][:,:,:]
    #nsr
    if "nsr_cu" in f:
        del f["nsr_cu"]
    if "nsr_cu_average" in f:
        del f["nsr_cu_average"]

    dset_nsr_cu = f.create_dataset("nsr_cu", (n_geomet,n_simulations,nblocks), dtype='f')
    dset_average_nsr_cu = f.create_dataset("nsr_cu_average", (nblocks,), dtype='f')

    #Au
    #grades
    dset_grades_au = f["simulation_au"][:,:]
    #concentrates
    dset_concentrates_au = f["geomet_au_con"][:,:,:]
    #recovery
    dset_recovery_au = f["geomet_au_rec"][:,:,:]
    #nsr
    if "nsr_au" in f:
        del f["nsr_au"]
    if "nsr_au_average" in f:
        del f["nsr_au_average"]
    dset_nsr_au = f.create_dataset("nsr_au", (n_geomet,n_simulations,nblocks), dtype='f')
    dset_average_nsr_au = f.create_dataset("nsr_au_average", (nblocks,), dtype='f')

    #F
    #concentrates

    #NSR total
    if "nsr" in f:
        del f["nsr"]
    if "nsr_average" in f:
        del f["nsr_average"]
    dset_nsr = f.create_dataset("nsr", (n_geomet,n_simulations,nblocks), dtype='f')
    dset_average_nsr = f.create_dataset("nsr_average", (nblocks,), dtype='f')


    grade = np.empty((n_simulations,2))
    concentrate = np.empty((n_geomet,2))
    recovery = np.empty((n_geomet,2))

    dset_nsr_cu_tmp = np.empty((n_geomet,n_simulations))
    dset_nsr_au_tmp = np.empty((n_geomet,n_simulations))

    grades_tmp = np.empty((nblocks*n_simulations,2))
    concentrate_tmp = np.empty((nblocks * n_simulations,2))
    recovery_tmp = np.empty((nblocks * n_simulations,2))

    for k in xrange(n_geomet):
        grades_cu_tmp = dset_grades_cu[:,:]
        grades_au_tmp = dset_grades_au[:,:]
        concentrate_cu_tmp = dset_concentrates_cu[k,:,:]
        concentrate_au_tmp = dset_concentrates_au[k,:,:]
        recovery_cu_tmp = dset_recovery_cu[k,:,:]
        recovery_au_tmp = dset_recovery_au[k,:,:]
            
        grade_cu = grades_cu_tmp.flatten()
        grade_au = grades_au_tmp.flatten()
        recovery_cu = recovery_cu_tmp.flatten()
        recovery_au = recovery_au_tmp.flatten()
        concentrate_cu = concentrate_cu_tmp.flatten()
        concentrate_au = concentrate_au_tmp.flatten()

        metal_node = root["metals"][0]
        price_cu = metal_node["price"]
        price_to_ton_cu = metal_node["unit2ton"]
        deduct_cu = metal_node["deduction"]
        deduct_proportion_cu = metal_node["payableRate"]
        refining_cost_cu = metal_node["refiningCharge"] #US$/lb refiningCharge

        metal_node = root["metals"][1]
        price_au = metal_node["price"]
        price_to_ton_au = metal_node["unit2ton"]
        deduct_au = metal_node["deduction"]
        deduct_proportion_au = metal_node["payableRate"]
        refining_cost_au = metal_node["refiningCharge"] #US$/oz

        smelter_charges = root["smelterCharge"]
        freight_charges = root["freightCharge"]
        penalties = root["penalties"]
        verbose = False
        #marketing_charge = root["marketing"]

        ret = nsr.calculateNSRLane(
        #cu
        grade_cu,
        recovery_cu,
        concentrate_cu,
        price_cu,
        price_to_ton_cu,
        #au
        grade_au,
        recovery_au,
        concentrate_au,
        price_au,
        price_to_ton_au,
        #cost cu
        deduct_cu,
        deduct_proportion_cu,
        refining_cost_cu,
        #cost au
        deduct_au,
        deduct_proportion_au,
        refining_cost_au,
        #common costs
        smelter_charges,
        freight_charges,
        verbose = False)        
        
        nsr_ore_cu,nsr_concentrate_cu,nsr_ore_au,nsr_concentrate_au = ret
        
        print nsr_ore_cu.shape,nsr_concentrate_cu.shape,nsr_ore_au.shape,nsr_concentrate_au.shape
        
        print grade_cu[:10]
        print recovery_cu[:10]
        print concentrate_cu[:10]
        print nsr_ore_cu[:10]
        print nsr_concentrate_cu[:10]
        dmt = grade_cu * (recovery_cu/100.0) / concentrate_cu
        print dmt[:10]
        
        #nsr_cu = nsrConcentrate[:,0].reshape(n_simulations,nblocks)
        #nsr_au = nsrConcentrate[:,1].reshape(n_simulations,nblocks)

        nsr_cu = nsr_ore_cu.reshape(n_simulations,nblocks)
        nsr_au = nsr_ore_au.reshape(n_simulations,nblocks)

        dset_nsr_cu[k,:,:] = nsr_cu
        dset_nsr_au[k,:,:] = nsr_au

        dset_nsr[k,:,:] = nsr_cu + nsr_au

    dset_average_nsr_cu[:] = np.mean(dset_nsr_cu[:,:,:],axis=(0,1))
    dset_average_nsr_au[:] = np.mean(dset_nsr_au[:,:,:],axis=(0,1))

    dset_average_nsr[:] = dset_average_nsr_cu[:] + dset_average_nsr_au[:]
