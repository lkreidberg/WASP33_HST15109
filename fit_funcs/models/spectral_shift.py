import sys
sys.path.insert(0,'..')
import numpy as np
from read_data import Data

def spectral_shift(t, data, params):
    xshift = params

    return 1. + xshift*data.spectral_shift
