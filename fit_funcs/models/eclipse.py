import batman

def eclipse(t, p, data, v_num):
    params = batman.TransitParams()

    params.t0 = p[data.par_order['t0']*data.nvisit + v_num]
    params.t_secondary = p[data.par_order['t_secondary']*data.nvisit + v_num]
    params.per = p[data.par_order['per']*data.nvisit + v_num]
    params.rp = p[data.par_order['rp']*data.nvisit + v_num]
    params.a = p[data.par_order['a']*data.nvisit + v_num]
    params.inc = p[data.par_order['inc']*data.nvisit + v_num]
    params.ecc = p[data.par_order['ecc']*data.nvisit + v_num]
    params.w = p[data.par_order['w']*data.nvisit + v_num]
    params.u = [p[data.par_order['u1']*data.nvisit + v_num], p[data.par_order['u2']*data.nvisit + v_num]]
    params.fp = p[data.par_order['fp']*data.nvisit + v_num]
    params.limb_dark = "quadratic"		#FIXME - specify this value in one of the config files
    
    m = batman.TransitModel(params, t, supersample_factor=3, exp_time = data.exp_time/24./60./60., transittype="secondary")
    return m.light_curve(params)
