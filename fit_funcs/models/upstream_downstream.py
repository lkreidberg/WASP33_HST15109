import sys
sys.path.insert(0,'..')
import numpy as np
from read_data import Data

def upstream_downstream(idx, data, params):
    scale = params
    return 1. + scale*data.scan_direction[idx]
