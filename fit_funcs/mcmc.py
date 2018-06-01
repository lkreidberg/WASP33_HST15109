import emcee

def format_params_for_mcmc(params, obs_par, fit_par):	#FIXME: make sure this works for cases when nvisit>1
    nvisit = int(obs_par['nvisit'])				
    theta = []

    for i in range(len(fit_par)):
            if fit_par['fixed'][i].lower() == "false":
                    if fit_par['tied'][i].lower() == "true": theta.append(params[i*nvisit])
                    else: 
                            for j in range(nvisit): theta.append(params[i*nvisit+j])
    return np.array(theta)

def mcmc_output(samples, params, obs_par, fit_par, data):	#FIXME: make sure this works for cases when nvisit>1
    nvisit = int(obs_par['nvisit'])				
    labels = []

    for i in range(len(fit_par)):
            if fit_par['fixed'][i].lower() == "false":
                    if fit_par['tied'][i].lower() == "true": labels.append(fit_par['parameter'][i])
                    else: 
                            for j in range(nvisit): labels.append(fit_par['parameter'][i]+str(j))
    fig = corner.corner(samples, labels=labels, show_titles=True)
    current_time = datetime.now().time()
    figname = "pairs_"+current_time.isoformat()+".png"
    fig.savefig(figname)


def format_params_for_Model(theta, params, obs_par, fit_par):
    nvisit = int(obs_par['nvisit'])
    params_updated = []
    iter = 0									#this should be a more informative name FIXME
    for i in range(len(fit_par)):
            if fit_par['fixed'][i].lower() == "true": 
                    for j in range(nvisit): 
                            params_updated.append(params[i*nvisit+j])
            else:
                    if fit_par['tied'][i].lower() == "true": 
                            for j in range(nvisit): params_updated.append(theta[iter])
                            iter += 1
                    else: 
                            for j in range(nvisit): 		
                                    params_updated.append(theta[iter])
                                    iter += 1
    return np.array(params_updated)

def mcmc_fit(file_name, obs_par, fit_par, flags):
    data = Data(file_name, obs_par, fit_par)
    if flags['divide-white']:
            sys_vector = np.genfromtxt("white_systematics.txt")
            data.all_sys = sys_vector
            data.nfree_param -= 2
            data.dof += 2
            print "subtracting 2 from dof for divide-white"


    print "e1", np.median(data.err)
    data, params = least_sq_fit(file_name, obs_par, fit_par, data, flags, myfuncs)		#starting guess
    print "e2", np.median(data.err)
    theta = format_params_for_mcmc(params, obs_par, fit_par)	

    ndim, nwalkers = len(theta), 50					#FIXME set nwalkers is a config file
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args = (params, data, obs_par, fit_par, flags))

    pos = [theta + 1e-5*np.random.randn(ndim) for i in range(nwalkers)]

    sampler.run_mcmc(pos,1000)
    #sampler.run_mcmc(pos,10000)
    pickle.dump([data, params, sampler.chain], open("mcmc_out."+"{0:0.2f}".format(data.wavelength)+".p", "wb"))

    samples = sampler.chain[:, 100:, :].reshape((-1, ndim))
    mcmc_output(samples, params, obs_par, fit_par, data)

    medians = []
    errors = []

    for i in range(len(theta)):
            q = quantile(samples[:, i], [0.16, 0.5, 0.84])
            medians.append(q[1])
            errors.append(q[2] - q[1])
    return data.wavelength, medians[0], errors[0], samples


def lnprior(theta, params, data, obs_par, fit_par, flags):
    #updated_params = format_params_for_Model(theta, params, obs_par, fit_par)
    #if sum(np.array(updated_params[data.par_order['c']*data.nvisit:(1 + data.par_order['c'])*data.nvisit])<0)>0: return -np.inf
    #else: return 0.
    return 0.
    

def lnprob(theta, params, data, obs_par, fit_par, flags):
    updated_params = format_params_for_Model(theta, params, obs_par, fit_par)
    #plot(updated_params, data, flags)
    m = Model(updated_params, data, flags, myfuncs)
    lp = lnprior(theta, params, data, obs_par, fit_par, flags, myfuncs)
    #cs = updated_params[data.par_order['c']*data.nvisit:(1 + data.par_order['c'])*data.nvisit]
    #print cs[0], cs[1], m.ln_likelihood
    return m.ln_likelihood + lp
