import sys

sys.path += ["/home/esepulveda/Documents/projects/geostatpy"]


import numpy as np
import csv
import os.path

if __name__ == "__main__":
    model_block_size = 80 * 70 * 100
    nvar = 6
    nsim = 50
    results_path = "/home/esepulveda/Documents/projects/newcrest/optimisation/sgsim/results-2"
    
    simulation_results = ["sgsim-v{var}.out".format(var=i+1) for i in xrange(nvar)]

    data = np.empty((model_block_size * nsim,nvar))
    
    for i,sim in enumerate(simulation_results):
        print "processing", sim
        with open(os.path.join(results_path,sim),"r") as fin:
            reader = csv.reader(fin,delimiter=" ",skipinitialspace = True)

            row = reader.next() #title
            row = reader.next() #variables and grid size
            check_m_block_size = int(row[1]) * int(row[2]) * int(row[3])
            row = reader.next() #variable name
            
            assert model_block_size == check_m_block_size
            
            for k,row in enumerate(reader):
                data[k,i] = float(row[0])
            

    header = '''Input data & PPMT transformed variables
6
PPMT:au_ppm
PPMT:cu_pct
PPMT:s_pct
PPMT:cucn_ppm
PPMT:fe_pct
PPMT:mo_ppm'''


    np.savetxt(os.path.join(results_path,"simulated_values.csv"),data,fmt="%.8f",delimiter=" ",header=header)

