import sys
sys.path.insert(0,'..')
import numpy as np
from read_data import Data

def spectral_shift(t, data, params):
    shift = params

    print 1. + shift*data.spectral_shift
    return 1. + shift*data.spectral_shift
