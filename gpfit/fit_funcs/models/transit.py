import batman

def transit(idx, data, params):
        t = data.time[idx]

	p = batman.TransitParams()
        t0, per, rp, a, inc, ecc, w, u1, u2, limb_dark = params
        #FIXME
        if limb_dark == 2: p.limb_dark = "quadratic"
        else: 
            print "unsupported limb darkening parameter"
            return 0

	p.t0 = t0
	p.per = per
	p.rp = rp
	p.a = a
	p.inc = inc
	p.ecc = ecc
	p.w = w
	p.u = [u1, u2]
	
	m = batman.TransitModel(
            p, t, supersample_factor=3, exp_time = data.exp_time/24./60./60.
        )
	return m.light_curve(p)
