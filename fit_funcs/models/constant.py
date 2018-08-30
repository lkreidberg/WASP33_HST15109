import sys
sys.path.insert(0,'..')
import numpy as np
from read_data import Data

def constant(idx, data, params):
    C = params
    return 1. + C*np.ones_like(data.time[idx])
