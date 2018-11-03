import sys
sys.path.insert(0,'..')
import numpy as np
from read_data import Data

def spatial_shift(t, data, params):
    yshift = params

    return 1. + yshift*data.spatial_shift
