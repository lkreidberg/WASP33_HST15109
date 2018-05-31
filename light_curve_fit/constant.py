import numpy as np

def constant(t, p, data, v_num):

    C = p[data.par_order['c']*data.nvisit + v_num]

    return 1.+C*np.ones_like(t)
