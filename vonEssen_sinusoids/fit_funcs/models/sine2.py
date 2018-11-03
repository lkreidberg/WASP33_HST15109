import sys
sys.path.insert(0,'..')
from read_data import Data
import numpy as np

def sine2(idx, data, params):
    a1, omega1, phi1, a2, omega2, phi2, a3, omega3, phi3, a4, phi4, omega4 = params
#    print a1, a2, a3, omega1, omega2, omega3
    TWOPI = np.pi*2.

    return  ( 1. + a1*np.sin(TWOPI*omega1*data.t_vis[idx] + phi1)  +
              + a2*np.sin(TWOPI*omega2*data.t_vis[idx] + phi2) 
              + a3*np.sin(TWOPI*omega3*data.t_vis[idx] + phi3) 
              + a4*np.sin(TWOPI*omega4*data.t_vis[idx] + phi4) 
         )


    
