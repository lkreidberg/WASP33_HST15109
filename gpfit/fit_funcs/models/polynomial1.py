import sys
sys.path.insert(0,'..')
from read_data import Data
import numpy as np

def polynomial1(idx, data, params):
    v = params
    #FIXME data.t_vis won't work if there are multiple visits
    return 1. + v*data.t_vis[idx]
