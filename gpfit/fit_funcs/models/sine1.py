import sys
sys.path.insert(0,'..')
from read_data import Data
import numpy as np

def sine1(idx, data, params):
    a1, omega1, phi1  = params
    return 1. + a1*np.sin(omega1*data.t_vis[idx] + phi1)
