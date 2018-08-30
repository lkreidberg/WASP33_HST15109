import sys
sys.path.insert(0, './fit_funcs/')
sys.path.insert(0, 'fit_funcs/models')
import pickle
import os 
import glob
import numpy as np
from read_data import Data
from model import Model
import matplotlib.pyplot as plt

def quantile(x, q): return np.percentile(x, [100. * qi for qi in q]) 

#path = "mcmc_output_ackbar_25bins/"
path = "mcmc_output_mr_25bins/"

mcmc = glob.glob(os.path.join(path, "mcmc*.p"))
lsq = glob.glob(os.path.join(path, "lsq*.p"))

print "FIXME setting ndim by hand"
print "AHHHHHHHHHHH"

if path == "mcmc_output_ackbar_25bins/": 
    ndim = 10
elif path == "mcmc_output_mr_25bins/": 
    ndim = 9
else: ndim = 0

f = open(path + "emission_spectrum.txt", "w")

xs, ys, es = [], [], []

for m, l in zip(mcmc, lsq):
    data, params, chain = pickle.load(open(m, "r"))
    data, model = pickle.load(open(l, "r"))


    samples = chain[:, 500:, :].reshape((-1, ndim))
    
    medians, errors = [], []

    for i in range(ndim):
            q = quantile(samples[:, i], [0.16, 0.5, 0.84])
            medians.append(q[1])
            errors.append(q[2] - q[1])
    
    y, ye = medians[0], errors[0]
    print data.wavelength,  y, ye 
    print>>f, data.wavelength,  y, ye 

    xs.append(data.wavelength)
    ys.append(y)
    es.append(ye)

f.close()

plt.errorbar(np.array(xs), np.array(ys), yerr = np.array(es), fmt = '.k')
plt.show()
