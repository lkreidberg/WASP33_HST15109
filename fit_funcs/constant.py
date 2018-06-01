import numpy as np

def constant(t, params):
    C = params
    return 1.+C*np.ones_like(t)
