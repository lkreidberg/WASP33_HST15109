import sys
sys.path.insert(0,'..')
from read_data import Data
import numpy as np

def essensine(idx, data, params):
    a1, phi1, a2, phi2, a3, phi3 = params

    omega1, omega2, omega3 = 20.1621, 21.0606, 9.8436
    TWOPI = 2.*np.pi

    return  (1. + a1*np.sin(TWOPI*omega1*data.t_vis[idx] + phi1) +  
                  a2*np.sin(TWOPI*omega2*data.t_vis[idx] + phi2) +
                  a3*np.sin(TWOPI*omega3*data.t_vis[idx] + phi3)
            ) 


    
