import sys
sys.path.insert(0,'..')
from read_data import Data
import numpy as np

def model_ramp(t, data, params):
    r1, r2, r3 = params
    #FIXME data.t_vis won't work if there are multiple visits
    return 1.0 - np.exp(-r1*data.t_orb - r2 - r3*data.t_delay)
