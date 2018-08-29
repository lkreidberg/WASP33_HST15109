import sys
sys.path.insert(0, './fit_funcs/')
sys.path.insert(0, 'fit_funcs/models')
import pickle
import os 
import glob
import numpy as np
from read_data import Data
from model import Model

def quantile(x, q): return np.percentile(x, [100. * qi for qi in q]) 

path = "mcmc_output/21bins/"

mcmc = glob.glob(os.path.join(path, "mcmc*.p"))
lsq = glob.glob(os.path.join(path, "lsq*.p"))

print "FIXME setting ndim by hand"

print "#wavelength, depth, depth_err, u1, u1_err, chi2red"
for m, l in zip(mcmc, lsq):
    data, params, chain = pickle.load(open(m, "r"))
    data, model = pickle.load(open(l, "r"))

    ndim = 5
    samples = chain[:, 500:, :].reshape((-1, ndim))
    
    medians, errors = [], []

    for i in range(ndim):
            q = quantile(samples[:, i], [0.16, 0.5, 0.84])
            medians.append(q[1])
            errors.append(q[2] - q[1])
    
    edepth, edepth_err =  (medians[0])**2, 2.*errors[0]*medians[0] 
    u1, u1_err =  medians[1], errors[1]
    print data.wavelength, "{0:0.6f}".format(depth), \
            "{0:0.1e}".format(depth_err*max(1., np.sqrt(model.chi2red))), \
            "{0:0.3f}".format(u1), "{0:0.3f}".format(u1_err), \
            "{0:0.2f}".format(model.chi2red)
