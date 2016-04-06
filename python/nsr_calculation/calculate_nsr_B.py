'''This case is considering only grade uncertainty. 50 simulation of grades and one best geomet model
'''
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

nsr_path = "/home/esepulveda/Documents/projects/newcrest/optimisation/nsr_results"
base_path = "/home/esepulveda/Documents/projects/newcrest/optimisation/results"

reload_density = False

#define or load grid
blockmodel3D = None
if os.path.exists(os.path.join(base_path,"grid3d.dump")):
    with open(os.path.join(base_path,"grid3d.dump"),"r") as fin: 
        try:
            blockmodel3D = pickle.load(fin)
            print "block model definition loaded"
        except Exception as e:
            print "block model definition loaded",e
            blockmodel3D = None

nx,ny,nz = blockmodel3D.nodes

nsim = 50
nblocks = nx*ny*nz

#create nsr dataset
nsr_filename = os.path.join(nsr_path,"nsr_B.mat")
h5nsr = h5py.File(nsr_filename, "a")

if True:
    dset_nsr = h5nsr.create_dataset("nsr", (nblocks,n_simulations), dtype=np.float32)
else:
    dset_nsr = h5nsr["nsr"]

home = "/home/esepulveda/Documents/projects/newcrest/optimisation"

json_parameters = os.path.join(home,"data","pc1s1.json");
root = json.load(open(json_parameters))

#grades
h5grades = h5py.File(os.path.join(home,"results","grade_simulations.mat"),"r")

#grades datasets
dset_grades_cu = h5grades["simulations/cu"]
dset_grades_au = h5grades["simulations/au"]

#geomet datasets
h5geomet = h5py.File(os.path.join(home,"geomet_results","geomet_B.mat"),"r")
dset_concentrates_cu = h5geomet["cu/concentrated"]
dset_recovery_cu = h5geomet["cu/recovery"]
dset_concentrates_au = h5geomet["au/concentrated"]
dset_recovery_au = h5geomet["au/recovery"]
dset_concentrates_f = h5geomet["f/concentrated"]

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

for k in xrange(nsim):
    grades_cu = dset_grades_cu[:,k]
    grades_au = dset_grades_au[:,k]
    concentrate_cu = dset_concentrates_cu[:,k]
    concentrate_au = dset_concentrates_au[:,k]
    recovery_cu = dset_recovery_cu[:,k]
    recovery_au = dset_recovery_au[:,k]

    nsr_ore_cu,nsr_concentrate_cu,nsr_ore_au,nsr_concentrate_au = nsr.calculateNSRLane(
        #cu
        grades_cu,
        recovery_cu,
        concentrate_cu,
        price_cu,
        price_to_ton_cu,
        #au
        grades_au,
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
    
    dset_nsr[:,k] = nsr_ore_cu + nsr_ore_au
