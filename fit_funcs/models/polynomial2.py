import sys
sys.path.insert(0,'..')
from read_data import Data
import numpy as np

def polynomial2(idx, data, params):
    v, v2 = params
    return 1. + v*data.t_vis[idx] + v2*(data.t_vis[idx])**2
