import sys
sys.path.insert(0, './light_curve_fit')
import mpfit
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
from astropy.io import ascii
import getopt
from datetime import datetime
import time as pythontime
import multiprocessing as mp
import emcee
import pickle
from scipy.stats import norm
from read_data import LightCurveData
from plot_data import plot_raw, plot_fit
from formatter import FormatParams, PrintParams
from light_curve_model import Model

def E(t, eta, f, Etot, tau, E0):
	return eta*f/(eta*f/Etot + 1./tau) + (E0 - eta*f/(eta*f/Etot + 1./tau))*np.exp(-(eta*f/Etot +1./tau)*t)


def computeRMS(data, maxnbins=None, binstep=1, isrmserr=False):
	# bin data into multiple bin sizes
	npts    = data.size
	if maxnbins == None:
		maxnbins = npts/10.
	binsz   = np.arange(1, maxnbins+binstep, step=binstep)
	nbins   = np.zeros(binsz.size)
	rms     = np.zeros(binsz.size)
	rmserr  = np.zeros(binsz.size)
	for i in range(binsz.size):
		nbins[i] = int(np.floor(data.size/binsz[i]))
		bindata   = np.zeros(nbins[i], dtype=float)
# bin data
# ADDED INTEGER CONVERSION, mh 01/21/12
		for j in range(int(nbins[i])):
			bindata[j] = data[j*binsz[i]:(j+1)*binsz[i]].mean()
		# get rms
		rms[i]    = np.sqrt(np.mean(bindata**2))
		rmserr[i] = rms[i]/np.sqrt(2.*int(nbins[i]))
	# expected for white noise (WINN 2008, PONT 2006)
	stderr = (data.std()/np.sqrt(binsz))*np.sqrt(nbins/(nbins - 1.))
	if isrmserr == True:
		return rms, stderr, binsz, rmserr
	else:
		return rms, stderr, binsz


def quantile(x, q):
        return np.percentile(x, [100. * qi for qi in q])



			
def usage():
	cmd = sys.argv[0]
	sys.stderr.write('Usage: python %s OPTION\n\n' % os.path.basename(cmd))
	sys.stderr.write(
		'Allowed OPTION flags are:\n'
		'  --show-plot      		displays fitted light curve plots\n'
		'  --run-mcmc      		runs MCMC starting from least-squares best fit parameters\n'
		'  --run-pb         		runs prayer bead analysis\n'
		'  --plot-raw-data		plots raw light curve separated by visit\n'   
		'  --plot-sys			plots light curve fit with visit-long systematics included\n'   
		'  --path PATHNAME		specifies PATHNAME to light curves to be fit (default = ./spec_lc/*)\n'
		'  --fit-white FILENAME		fits the white light curve stored in FILENAME (default = "lc_white.txt"\n'
		'  -v               		prints fit diagnostic information\n'
		'  -o               		saves fit output to file\n'
		'  -h               		lists instructions for usage\n'
		'\n')
	sys.exit(1)

def make_dict(table):
	return {x['parameter']: x['value'] for x in table}


def weighted_mean(data, err):            #calculates the weighted mean for data points data with std devs. err
	ind = err != 0.0
        weights = 1.0/err[ind]**2
        mu = np.sum(data[ind]*weights)/np.sum(weights)
        var = 1.0/np.sum(weights)
        return [mu, np.sqrt(var)]                


def residuals(params, data, flags, fjac=None):					
	return [0, Model(params, data, flags).resid/data.err]

def least_sq_fit(file_name, obs_par, fit_par, data, flags):
	nvisit = int(obs_par['nvisit'])
	npar = len(fit_par)*nvisit

	#initializes least squares fit parameters
	parinfo = [{'value':0, 'fixed':0, 'limited':[0,0,], 'limits':[0.0,0.0], 'step':0.0} for j in range(npar)]
	params_s = []

	for i in range(npar/nvisit):								#loops through params
		for j in range(nvisit):								#loops through visits
			parinfo[i*nvisit+j]['value'] = fit_par['value'][i]			#sets initial guess value	
			parinfo[i*nvisit+j]['step'] = 0.01*np.abs(fit_par['value'][i])		#sets parameter step size
			if i==1: parinfo[i*nvisit+j]['step'] = 0.00001				#makes small step for t0
			parinfo[i*nvisit+j]['fixed'] = fit_par['fixed'][i].lower() == "true"	#sets whether parameter varies
			if j>0 and fit_par['tied'][i].lower() == "true":
				parinfo[i*nvisit+j]['tied'] = 'p[{0}]'.format(nvisit*i)		#ties parameters to first visit value
			if fit_par['lo_lim'][i].lower() == "true": 				#puts lower limit on parameter
				parinfo[i*nvisit+j]['limited'][0] = True
				parinfo[i*nvisit+j]['limits'][0] = fit_par['lo_val'][i]
			if fit_par['hi_lim'][i].lower() == "true": 				#puts upper limit on parameter
				parinfo[i*nvisit+j]['limited'][1] = True
				parinfo[i*nvisit+j]['limits'][1] = fit_par['hi_val'][i]
			params_s.append(fit_par['value'][i])

	params_s = np.array(params_s)
	data = LightCurveData(file_name, obs_par, fit_par)
	
	if flags['plot-raw-data']: plot_raw(data)
	fa = {'data':data, 'flags':flags}

	if flags['divide-white']:
		sys_vector = np.genfromtxt("white_systematics.txt")
		data.all_sys = sys_vector
		data.nfree_param -= 2
		data.dof += 2
#		print "subtracting 2 from dof for divide-white"

	m = mpfit.mpfit(residuals, params_s, functkw=fa, parinfo = parinfo, quiet=True) 
	model = Model(m.params, data, flags)

	m = mpfit.mpfit(residuals, params_s, functkw=fa, parinfo = parinfo, quiet=True) 
	model = Model(m.params, data, flags)
	"""#rescale error bars based on chi2
	#print "rescaling error bars to get chi2red = 1"
	#print "scale factor = ", np.sqrt(model.chi2red)
	print data.wavelength, np.sqrt(model.chi2red)
	data.err = data.err*np.sqrt(model.chi2red)
	m = mpfit.mpfit(residuals, params_s, functkw=fa, parinfo = parinfo, quiet=True) 
	model = Model(m.params, data, flags)"""
	
	if flags['output']: 
		f = open(flags['out-name'], "a")
		#print>>f, "{0:0.3f}".format(data.wavelength), "{0:0.6f}".format(m.params[data.par_order['rp']*nvisit]), "{0:0.6f}".format(m.perror[data.par_order['rp']*nvisit]), "{0:0.2f}".format(Model(m.params, data, flags).chi2red), "{0:0.3f}".format(m.params[data.par_order['u1']*nvisit]), "{0:0.3f}".format(m.perror[data.par_order['u1']*nvisit])
		print>>f, "{0:0.3f}".format(data.wavelength), "{0:0.6f}".format(m.params[data.par_order['rp']*nvisit]), "{0:0.6f}".format(m.perror[data.par_order['rp']*nvisit]), "{0:0.3f}".format(m.params[data.par_order['u1']*nvisit]), "{0:0.3f}".format(m.params[data.par_order['u2']*nvisit]), "{0:0.2f}".format(Model(m.params, data, flags).chi2red)
		pickle.dump([data, model], open("white_lc_fit.p", "wb"))
		f.close()


	if flags['verbose']: 
		model = Model(m.params, data, flags)
		print "{0:0.3f}".format(data.wavelength), "{0:0.2f}".format(model.chi2red)#, m.params[data.par_order['amp1']*nvisit], m.perror[data.par_order['amp1']*nvisit]
		PrintParams(m, data)

	if flags['show-plot']: plot_fit(m.params, data, flags, obs_par, plot_sys=flags['plot-sys'])
	return  data, m.params

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
	data = LightCurveData(file_name, obs_par, fit_par)
	if flags['divide-white']:
		sys_vector = np.genfromtxt("white_systematics.txt")
		data.all_sys = sys_vector
		data.nfree_param -= 2
		data.dof += 2
		print "subtracting 2 from dof for divide-white"


	print "e1", np.median(data.err)
	data, params = least_sq_fit(file_name, obs_par, fit_par, data, flags)		#starting guess
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
	m = Model(updated_params, data, flags)
	lp = lnprior(theta, params, data, obs_par, fit_par, flags)
	#cs = updated_params[data.par_order['c']*data.nvisit:(1 + data.par_order['c'])*data.nvisit]
	#print cs[0], cs[1], m.ln_likelihood
	return m.ln_likelihood + lp

def main():
	#parses command line input
	try: opts, args = getopt.getopt(sys.argv[1:], "hov", ["help", "show-plot", "run-mcmc", "plot-raw-data", "plot-sys", "path=", "fit-white=", "divide-white"]) 
	except getopt.GetoptError: usage()

	#defaults for command line flags
	verbose, output, show_plot, run_mcmc, run_lsq, plot_raw_data, plot_sys, path, fit_white, divide_white = False, False, False, False, True, False, False, "spec_lc", False, False

	for o, a in opts:
		if o in ("-h", "--help"): usage()
		elif o == "-o": output = True
		elif o == "-v": verbose = True
		elif o == "--show-plot": show_plot = True
		elif o == "--run-mcmc": run_mcmc, run_lsq = True, False
		elif o == "--run-lsq": run_lsq = True
		elif o == "--plot-raw-data": plot_raw_data = True
		elif o == "--plot-sys": plot_sys = True
		elif o == "--path": path = a
		elif o == "--fit-white": fit_white, white_file = True, a
		elif o == "--divide-white": divide_white = True
		else: assert False, "unhandled option"

	flags = {'verbose': verbose, 'show-plot': show_plot, 'plot-raw-data': plot_raw_data, 'plot-sys': plot_sys, 'output': output, 'out-name': 'none.txt', 'run-lsq': run_lsq, 'run-mcmc': run_mcmc, 'divide-white': divide_white, 'fit-white': fit_white}		#put these in dictionary right away!!

	#reads in observation and fit parameters
	obs_par = make_dict(ascii.read("config/obs_par.txt", Reader=ascii.CommentedHeader))
	fit_par = ascii.read("config/fit_par.txt", Reader=ascii.CommentedHeader)		

	files = glob.glob(os.path.join(path, "*"))		
	if fit_white: files = glob.glob(white_file)

	flags['out-name'] = "fit_" + pythontime.strftime("%Y_%m_%d_%H:%M") + ".txt"

	for f in files:
		data = LightCurveData(f, obs_par, fit_par)
		if flags['divide-white']:
			sys_vector = np.genfromtxt("white_systematics.txt")
			data.all_sys = sys_vector
		data, params = least_sq_fit(f, obs_par, fit_par, data, flags)
		m = Model(params, data, flags)
	
		#outfile = open("white_systematics.txt", "w")
		#for i in range(len(m.all_sys)): print>>outfile, m.all_sys[i]
		#outfile.close()


		if flags['run-mcmc']:
			output = mcmc_fit(f, obs_par, fit_par, flags)
		"""for i in range(len(output)): 
			temp = "mcmc_out_"+'{0:0.2f}'.format(output[i][0])
			np.save(temp, output[i][3])
			print output[i][0], output[i][1], output[i][2]
		outfile = open("mcmc_output.txt", "a")
		for i in range(len(output)): print>>outfile, output[i][0], output[i][1], output[i][2]
		outfile.close()"""
				

if __name__ == '__main__':
	main()
