import sys

sys.path += ["/home/esepulveda/Documents/projects/fortran"]
sys.path += ["/home/esepulveda/Documents/projects/newcrest/scripts/python"]

import pickle
import numpy as np
from ppreg import ProjectionPursuitRegression
import scipy
import scipy.interpolate
import matplotlib.pyplot as plt
import csv
import statistics
import sklearn
import os.path

if __name__ == "__main__":
    
    np.random.seed(123456)
    
    home_path = "/home/esepulveda/Documents/projects/newcrest/"

    data = np.loadtxt(os.path.join(home_path,"databases/recovery/recovery_stages_9_10_11.csv"),delimiter=",",skiprows=1)

    au_ppm = data[:,1]
    cu_pct = data[:,2]
    mo_ppm = data[:,3]
    fe_pct = data[:,4]
    s_pct = data[:,5]
    cucn_pct = data[:,6]
    f_ppm = data[:,7]
    cu_s = data[:,8]
    mo_pct = data[:,3] / 10000

    au_recleaner = data[:,16]
    cu_recleaner = data[:,17]
    f_recleaner = data[:,20]

    au_rougher_recovery = data[:,26]
    cu_rougher_recovery = data[:,27]
    

    n = data.shape[0]

    inputs = [au_ppm,cu_pct,mo_pct,fe_pct,s_pct,cucn_pct,f_ppm]

    x = np.empty((len(inputs),n))
    for i,inp in enumerate(inputs):
        x[i,:] = inp[:]

    outputs = [au_recleaner,cu_recleaner,f_recleaner,au_rougher_recovery,cu_rougher_recovery]
    error_limit = [0.1,0.05,0.8,0.05,0.05]
    
    
    samples = 1000
    
    ret = np.empty((len(outputs),samples))
    
    models = []
    model_names = ["au_recleaner_concentrate","cu_recleaner_concentrate","f_recleaner_concentrate","au_rougher_recovery","cu_rougher_recovery"]

    #creating csv outputs
    output_filename = os.path.join(home_path,"optimisation/ppr_models_5/summary_ppr_model_{output}.csv")
    
    csv_writer = {}
    model_output = {}
    for output in model_names:
        csv_writer[output] = csv.writer(open(output_filename.format(output=output),"w"),delimiter=",")
        csv_writer[output].writerow(["bootstrap","r","r2","mean error","std error"])
        #model outputs
        model_output[output] = []
        
    for j in xrange(samples):
        indices = np.random.randint(n, size=n)
        indices = list(set(indices)) #keep uniques
        indices.sort()
        
        for i,output in enumerate(outputs):
            y = output.copy()

            model = ProjectionPursuitRegression(6,6)
            model.fit(x[:,indices],y[indices].reshape((1,len(indices))))

            prediction = model.predict(x)
            
            residuals = (y - prediction)
            mean_error = np.mean(residuals)
            std_error = np.std(residuals)
            
            r = np.corrcoef(output,prediction)[0,1]
            r2 = max(0,sklearn.metrics.r2_score(y,prediction))
            
            row = [j,r,r2,mean_error,std_error]
            
            csv_writer[model_names[i]].writerow(row)
            
            #if np.abs(mean_error) <= error_limit[i]:
            model_output[model_names[i]] += [(abs(mean_error),r2,model,row)]

    #max_models = 10

    #for output in model_names:
    #    mo = model_output[output]
    #    assert len(mo) >= max_models,"bad max_models: " + output + ":" + str(len(mo))



    output_filename = os.path.join(home_path,"optimisation/ppr_models_7/summary_selected_ppr_models.csv")
    csv_writer = csv.writer(open(output_filename,"w"),delimiter=",")
    for output in model_names:
        mo = model_output[output]
        mo.sort()
        #mo.reverse()
        
        max_models = 50


        
        print output,len(mo)
        for i in xrange(max_models):
            abs_mean_error,r2,model,row = model_output[output][i]
            row.insert(0,output)
            row.insert(0,str(i+1))
            csv_writer.writerow(row)
            
            fn = os.path.join(home_path,"optimisation/ppr_models_7/ppr_model_{name}_{number}.dump")
            with open(fn.format(name=output,number=i),"w") as fout:
                pickle.dump(mo[i][2],fout,-1)


            print i,mo[i][0],mo[i][1]

