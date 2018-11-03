import sys
sys.path.insert(0,'..')
from read_data import Data
import numpy as np

def model_ramp(idx, data, params):
    r1, r2, r3  = params
    return 1.0 - np.exp(-r1*data.t_orb[idx] - r2 - r3*data.t_delay[idx])
