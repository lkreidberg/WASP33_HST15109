import sys
sys.path.insert(0,'..')
from read_data import Data
import numpy as np
#import matplotlib.pyplot as plt

def vonessen(idx, data, params):
    print params
    x, a1, a2, a3, a4, a5, a6, a7, a8, p1, p2, p3, p4, p5, p6, p7, p8, scale, q2, q3, q4, q5, q6, q7, q8 = params 
    

    omega1 = 20.1621
    omega2 = 21.0606
    omega3 = 8.8436
    omega4 = 24.8835
    omega5 = 20.5353
    omega6 = 34.1252
    omega7 = 8.3084
    omega8 = 10.8249

    a1 *= scale
    a2 *= scale
    a3 *= scale
    a4 *= scale
    a5 *= scale
    a6 *= scale
    a7 *= scale
    a8 *= scale

    TWOPI = np.pi*2.

    t_vis = data.t_vis[idx]

    return ( 1. + a1*np.sin(TWOPI*omega1*t_vis + p1) 
                + a2*np.sin(TWOPI*omega2*t_vis + p2) 
                + a3*np.sin(TWOPI*omega3*t_vis + p3) 
                + a4*np.sin(TWOPI*omega4*t_vis + p4) 
                + a5*np.sin(TWOPI*omega5*t_vis + p5) 
                + a6*np.sin(TWOPI*omega6*t_vis + p6) 
                + a7*np.sin(TWOPI*omega7*t_vis + p7) 
                + a8*np.sin(TWOPI*omega8*t_vis + p8) 
           )

       
