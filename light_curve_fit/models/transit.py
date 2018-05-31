import batman

def transit(t, params, data, v_num):
	p = batman.TransitParams()

        p.t0, p.t_secondary, p.per, p.rp, p.a, p.inc, p.ecc, p.w, p.u, p.limb_dark = params
        if p.limb_dark = 2: p.limb_dark = "quadratic"
        else: 
            print "unsupported limb darkening parameter"
            return 0

	"""bat_params.t0 = p.t0[v_num]
	bat_params.t_secondary = p.t_secondary[v_num]
	bat_params.per = p.per[v_num]
	bat_params.rp = p.rp[v_num]
	bat_params.a = p.a[v_num]
	bat_params.inc = p.inc[v_num]
	bat_params.ecc = p.ecc[v_num]
	bat_params.w = p.w[v_num]
	bat_params.u = [p.u1[v_num], p.u2[v_num]]
	bat_params.limb_dark = "quadratic"		#FIXME - specify this value in one of the config files"""
	
	m = batman.TransitModel(p, t, supersample_factor=3, exp_time = data.exp_time/24./60./60.)
	return m.light_curve(p)
