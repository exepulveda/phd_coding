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
nblocks = nx*ny*nz

#create nsr dataset
nsr_filename = os.path.join(ds_path,"nsr_average_inputs.mat")
h5nsr = h5py.File(nsr_filename, "r+")

if False:
    dset_nsr = h5nsr.create_dataset("average_input/nsr", (nblocks,n_simulations), dtype='f')
    dset_nsr_average = h5nsr.create_dataset("average_input/nsr_average", (nblocks,), dtype='f')

dset_nsr = h5nsr["average_input/nsr"]
dset_nsr_average = h5nsr["/average_input/nsr_average"]

home = "/home/esepulveda/Documents/projects/newcrest/optimisation"

json_parameters = os.path.join(home,"data","pc1s1.json");
root = json.load(open(json_parameters))

h5grades = h5py.File(os.path.join(home,"results","average_grade_simulations.mat"),"r")

#cu datasets
dset_grades_cu = h5grades["simulations/cu"]
dset_concentrates_cu = h5nsr["cu/concentrated"]
dset_recovery_cu = h5nsr["cu/recovery"]

#au datasets
dset_grades_au = h5grades["simulations/au"]
dset_concentrates_au = h5nsr["au/concentrated"]
dset_recovery_au = h5nsr["au/recovery"]

#fluorine
dset_concentrates_f = h5nsr["f/concentrated"]

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

    grade_cu = np.empty((e-s,n_simulations))
    grade_au = np.empty((e-s,n_simulations))
    
    grades_cu_tmp = dset_grades_cu[s:e]
    grades_au_tmp = dset_grades_au[s:e]
    concentrate_cu = dset_concentrates_cu[s:e,:]
    concentrate_au = dset_concentrates_au[s:e,:]
    recovery_cu = dset_recovery_cu[s:e,:]
    recovery_au = dset_recovery_au[s:e,:]
    
    #broadcast grades
    for i in xrange(n_simulations):
        grade_cu[:,i] = grades_cu_tmp[:]
        grade_au[:,i] = grades_au_tmp[:]

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
    
    
