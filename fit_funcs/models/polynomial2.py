import sys
sys.path.insert(0,'..')
from read_data import Data
import numpy as np

def polynomial2(t, data, params):
    v, v2 = params
    #FIXME data.t_vis won't work if there are multiple visits
    return 1. + v*data.t_vis + v2*(data.t_vis)**2
