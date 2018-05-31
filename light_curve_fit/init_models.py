import numpy as np

models = []
param_indices = []

def init_models(params, data, models_from_run):
    #self.per = params[data.par_order['per']*data.nvisit:(1 + data.par_order['per'])*data.nvisit]

    po = data.par_order

    for m in models_from_run:
        models += m

        if m == "constant":
            param_indices += [po['C']*data.nvisit]

        #if m == "transit":
        #    param_indices += [po['t0'], po['t_secondary'], po['per'], po['rp'], 
        #           po['a'], po['inc'], po['ecc'], po['w'], po['u'], 
        #          po['limb_dark']]
            

    
            
