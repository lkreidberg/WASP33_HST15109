import sys
sys.path.insert(0,'..')
from read_data import Data
import numpy as np

def sine2(idx, data, params):
    a1, omega1, phi1, a2, omega2, phi2 = params
    lc = ( 1. + a1*np.sin(omega1*data.t_vis[idx] + phi1) ) * \
           ( 1. + a2*np.sin(omega2*data.t_vis[idx] + phi2) )

    #np.save("bestfit_sine", [data.time, lc])
    return ( 1. + a1*np.sin(omega1*data.t_vis[idx] + phi1) ) * \
           ( 1. + a2*np.sin(omega2*data.t_vis[idx] + phi2) )

    
