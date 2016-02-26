import sys

sys.path += ["/home/esepulveda/Documents/projects/fortran"]

import numpy as np
import sys
import h5py
import ppr_model
import pickle
import os.path

def create_distribution(mean_values, dimens, minvalue = 0.0, maxvalue = np.inf):
    print "create_distribution..."
    values = np.random.normal(size=dimens)
    stdcu = np.std(mean_values)
    
    print "stdcu",stdcu
    
    #rescaling
    values = values * (stdcu)  + mean_values

    #clipping
    values = np.clip(values,minvalue,maxvalue)

    print "create_distribution...DONE"
    
    return values
    

class DrawPoint(object):
    def __init__(self,x,y):
        self._offset = [x,y]
        
    def stream(self):
        n_windows = 4

        i_j = []
        for i in xrange(-size_x,size_x):
            for j in xrange(-size_y,size_y):
                i_j = (i,j)
        
        ret = []
        for k in xrange(self.hight):
            for i,j in i_j:
                ret += [(i,j,k)]
                
        return ret

def recovery(x,a,b):
    #return np.min(np.max(a*np.log(x+0.001) + b,0.99),0.0)
    ret = a*np.log(x+0.5) + b

    indices = np.where(ret>= 1.0)[0]
    
    #print len(indices),ret[indices]
    #if len(indices) > 10:
    #    quit()
    ret[indices] = 0.99
    
    return ret




nx = 15
ny = 30
nz = 50
nblocks = nx*ny*nz
n_simulations = 25
n_geomet = 25

#load block_model example
blockmodel_data = np.loadtxt("../data/block_model_22500.csv",delimiter=",",skiprows=1)

#centroid_x,centroid_y,centroid_z,dim_x,dim_y,dim_z,volume,au_ok,cu_ok,mo_ok,f_ok,bulk_density,s_ok,ag_ok,fe_ok,f_con_avg,cu_con_grade,rec_cu,rec_au,grav_au_rec

n = len(blockmodel_data)

au_grades = blockmodel_data[:,7]
cu_grades = blockmodel_data[:,8]
mo_grades = blockmodel_data[:,9]
fe_grades = blockmodel_data[:,14]
s_grades = blockmodel_data[:,12]
f_grades = blockmodel_data[:,10]

density = blockmodel_data[:,11]
cu_concentrate = blockmodel_data[:,16]
cu_recovery = blockmodel_data[:,17] 
au_recovery = blockmodel_data[:,18]

#input data [au_ppm,cu_pct,mo_ppm,fe_pct,s_pct,f_ppm]
n_inputs = 6
input_data = np.empty((n,n_inputs))
input_data[:,0] = au_grades
input_data[:,1] = cu_grades
input_data[:,2] = mo_grades
input_data[:,3] = fe_grades
input_data[:,4] = s_grades
input_data[:,5] = f_grades

#minvalues = np.array([0.01,0.01,5.0,0.01,0.01,0.01])
minvalues = np.min(input_data,axis=0)
maxvalues = np.max(input_data,axis=0)*1.2
assert (minvalues.shape[0] == n_inputs)
assert (maxvalues.shape[0] == n_inputs)


#create covariance matrix
cov_matrix = np.cov(input_data,rowvar=0)
means = np.mean(input_data,axis=0)

print "...creating hdf5 datasets..."

#create hdf5 datasets
input_names = ["au","cu","mo","fe","s","f"]
input_variables = { "au": au_grades,"cu": cu_grades,"mo": mo_grades,"fe": fe_grades,"s": s_grades,"f": f_grades,}

with h5py.File("../data/bm_{0}_{1}_{2}_dataset.hdf5".format(nblocks,n_simulations,n_geomet), "w") as f:
    for k,v in input_variables.iteritems():
        dset = f.create_dataset("input_{}".format(k), (nblocks,), dtype='f')
        dset[:] = v

        f.create_dataset("simulation_{}".format(k), (n_simulations,nblocks), dtype='f')

    output_names = ["au_con","cu_con","f_con","au_rec","cu_rec"]
    for v in output_names:
        # geomet models
        f.create_dataset("geomet_{}".format(v), (n_geomet,n_simulations,nblocks), dtype='f')

    #density
    dset = f.create_dataset("density", (nblocks,), dtype='f')
    dset[:] = density

    #sample simulations
    print "...creating simulations..."

    simulation_samples = np.empty((n_inputs,n_simulations,nblocks))

    for i in xrange(n):
        inputs = input_data[i,:]
        samples = np.random.multivariate_normal(inputs,cov_matrix,size=(n_simulations,))

        for j in xrange(n_inputs):
            #clip
            simulation_samples[j,:,i] = np.clip(samples[:,j],minvalues[j],maxvalues[j])

    for j in xrange(n_inputs):
        #save simulations
        f["simulation_{}".format(input_names[j])][:,:] = simulation_samples[j,:,:]

    print "...running geomet models..."

    #run ppr to au_recovery
    model_names = ["au_recleaner_concentrate","cu_recleaner_concentrate","f_recleaner_concentrate","au_recleaner_recovery","cu_recleaner_recovery"]
    minimum_models = np.array([3.0,10.0,100.0,30.0,40.0])
    n_outputs = len(model_names)
    path = "/data/Linux/Documents/projects/newcrest/scripts/python"

    #load models
    ppr_models = {}
    for nm in model_names:
        if nm not in ppr_models:
            ppr_models[nm] = []
        for i in xrange(n_geomet):
            with open(os.path.join(path,"ppr_model_{}_{}.dump".format(nm,i)),"r") as fin:
                ppr_models[nm] += [pickle.load(fin)]

    geomet_simulations = np.empty((n_outputs,n_geomet,n_simulations,nblocks))

    inputs_tmp = np.empty((n_inputs,n_simulations*nblocks)) #,n_simulations,nblocks

    for j in xrange(n_geomet):
        for k,nm in enumerate(model_names):
            
            for i in xrange(n_inputs):
                inputs_tmp[i,:] = simulation_samples[i,:,:].flatten()

            georet = ppr_models[nm][j].predict(inputs_tmp)
            #clip
            georet = np.clip(georet,minimum_models[k],np.inf)
            geomet_simulations[k,j,:,:] = georet.reshape(n_simulations,nblocks)

    print "...saving geomet models..."

    for j in xrange(n_outputs):
        #print geomet_simulations[j,:,:,:].shape
        f["geomet_{}".format(output_names[j])][:,:,:] = geomet_simulations[j,:,:,:]

    print "...DONE..."

    quit()

    for i in xrange(nblocks):
        if (i % 100) == 0:
            print "calculing geomet to",i
        
        inputs = simulation_samples[:,:,i]
        #print inputs.shape

        for j in xrange(n_geomet):
            for k,nm in enumerate(model_names):
                #print ppr_models[nm][j].predict(inputs).shape
                georet = ppr_models[nm][j].predict(inputs)
                #clip
                georet = np.clip(georet,0,np.inf)
                geomet_simulations[k,j,:,i] = georet

            #geomet_simulations[k,:,j,i] = result

    print "...saving geomet models..."

    for j in xrange(n_outputs):
        #print geomet_simulations[j,:,:,:].shape
        f["geomet_{}".format(output_names[j])][:,:,:] = geomet_simulations[j,:,:,:]

    print "...DONE..."
