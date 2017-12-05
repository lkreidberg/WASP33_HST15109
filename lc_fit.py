import sys
#sys.path.insert(1, '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/')
import mpfit
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import glob
import os
from astropy.io import ascii
import getopt
import seaborn as sns
import matplotlib
from datetime import datetime
import time as pythontime
import multiprocessing as mp
import emcee
import batman
import pickle
#import corner
from scipy.stats import norm

sns.set_context("talk")
sns.set_style("white")
sns.set_style("ticks", {"xtick.direction":"in", "ytick.direction":"in"})
matplotlib.rcParams.update({'lines.markeredgewidth':0.3})
matplotlib.rcParams.update({'axes.formatter.useoffset':False})

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

def parse_skip_orbs(x):
	n = len(x)
	i = 0
	temp = 0
	skip_orbs = []
	while(i < n):
		if x[i] == "_": 
			skip_orbs.append(int(x[temp:i]))
			temp = i+1
		i += 1
	return np.array(skip_orbs)

class LightCurveData:
	"""
	doc
	"""
	def __init__(self, data_file, obs_par, fit_par):
		par_order = {line['parameter']: i for i, line in enumerate(fit_par)}

		d = np.genfromtxt(data_file)

		d = d[np.argsort(d[:,5])]				#sorts data by time

		use_first_exp = False
		if use_first_exp == False:
			ind = np.diff(d[:,5]) < 30./60./24.		#removes first exposure from each orbit
			d = d[1:][ind]

		orb_num = np.zeros_like(d[:,7])					#removes first orbit from each visit 
		orb = 0
		for i in range(1, len(orb_num)):
			if (d[i,5] - d[i-1,5]) > 0.5/24.: orb += 1	#stores data as separate visit if gap is longer than 9 hrs; FIXME: what if visits are closer together than 9 hours?
			orb_num[i] = orb


                #LK change - keep zeroth orbit! 12/1/17
		"""ind = orb_num == 0
		d = d[~ind]
		orb_num = orb_num[~ind]
		orb_num -= 1"""

		if obs_par['lc_type'] == "transit": skip_orbs = parse_skip_orbs(obs_par['skip_orbs_transit'])
		elif obs_par['lc_type'] == "eclipse": skip_orbs = parse_skip_orbs(obs_par['skip_orbs_eclipse'])
		elif obs_par['lc_type'] == "phase_curve": skip_orbs = parse_skip_orbs(obs_par['skip_orbs_phase_curve'])
		else: raise Exception("Unsupported lc_type")

		
		n = len(d)
		vis_num = np.zeros(n)
		t_vis = np.zeros(n) 
		t_orb = np.zeros(n) 
		t_delay = np.zeros(n) 
		
		visit = 0
		for i in range(1,n):
			vis_num[i - 1] = visit
			if (d[i,5] - d[i-1,5]) > 9./24.: visit += 1	#stores data as separate visit if gap is longer than 9 hrs; FIXME: what if visits are closer together than 9 hours?
		vis_num[-1] = visit

		nvisit = int(obs_par['nvisit'])
		norbit = int(obs_par['norb'])-1

		for i in range(nvisit):
			ind = vis_num == i
			t_vis[ind] = d[ind,5] - d[ind,5][0]
		
		orbs_per_visit = norbit/nvisit

		for i in range(norbit):
			ind = orb_num == i
			t_orb[ind] = d[ind,5] - d[ind,5][0]
			if i%orbs_per_visit == False: t_delay[ind] = 1.
	

		err = np.sqrt(d[:,2])
		flux = d[:,1]
		time  = d[:,5]
		scan_direction = d[:,8]

		wavelength = d[0,3]
		#wavelength = 0.7
		#print "setting wavelength by hand to fix_ld for white lc"

		#fixes limb darkening if "fix_ld" parameter is set to True in obs_par.txt
		if obs_par['fix_ld'].lower() == "true":
			ld = np.genfromtxt(obs_par['ld_file'])
			
			i = 0
			while(wavelength > ld[i,1]): i += 1
			
			u1 = ld[i, 3]
			u2 = ld[i, 4] 
			fit_par['value'][np.where(fit_par['parameter']=='u1')] = u1
			fit_par['fixed'][np.where(fit_par['parameter']=='u1')] = "true"
			fit_par['value'][np.where(fit_par['parameter']=='u2')] = u2
			fit_par['fixed'][np.where(fit_par['parameter']=='u2')] = "true"

		nfree_param = 0
		for i in range(len(fit_par)):
			if fit_par['fixed'][i].lower() == "false":
				if fit_par['tied'][i].lower() == "true": nfree_param += 1
				else: nfree_param += nvisit

		self.time = time
		self.flux = flux
		self.err = err
		self.wavelength = wavelength
		self.exp_time = float(obs_par['exp_time'])
		self.nvisit = nvisit
		self.vis_num = vis_num
		self.orb_num = orb_num
		self.scan_direction = scan_direction
		self.t_vis = t_vis
		self.t_orb = t_orb
		self.t_delay = t_delay
		self.par_order = par_order
		self.nfree_param = nfree_param
		self.dof = n - nfree_param 
		self.npoints = n
		self.lc_type = obs_par['lc_type']
		self.all_sys = None
		self.u1 = 0.
		self.u2 = 0.


class FormatParams: 
	"""
	doc
	"""
	def __init__(self, params, data):
		self.per = params[data.par_order['per']*data.nvisit:(1 + data.par_order['per'])*data.nvisit]
		self.t0 = params[data.par_order['t0']*data.nvisit:(1 + data.par_order['t0'])*data.nvisit]
		self.t_secondary = params[data.par_order['t_secondary']*data.nvisit:(1 + data.par_order['t_secondary'])*data.nvisit]
		self.w = params[data.par_order['w']*data.nvisit:(1 + data.par_order['w'])*data.nvisit]
		self.a = params[data.par_order['a']*data.nvisit:(1 + data.par_order['a'])*data.nvisit]
		self.inc = params[data.par_order['inc']*data.nvisit:(1 + data.par_order['inc'])*data.nvisit]
		self.rp = params[data.par_order['rp']*data.nvisit:(1 + data.par_order['rp'])*data.nvisit]
		self.fp = params[data.par_order['fp']*data.nvisit:(1 + data.par_order['fp'])*data.nvisit]
		self.u1 = params[data.par_order['u1']*data.nvisit:(1 + data.par_order['u1'])*data.nvisit]
		self.u2 = params[data.par_order['u2']*data.nvisit:(1 + data.par_order['u2'])*data.nvisit]
		self.ecc = params[data.par_order['ecc']*data.nvisit:(1 + data.par_order['ecc'])*data.nvisit]
		self.c = params[data.par_order['c']*data.nvisit:(1 + data.par_order['c'])*data.nvisit]
		self.v = params[data.par_order['v']*data.nvisit:(1 + data.par_order['v'])*data.nvisit]
		self.v2 = params[data.par_order['v2']*data.nvisit:(1 + data.par_order['v2'])*data.nvisit]
		self.r1 = params[data.par_order['r1']*data.nvisit:(1 + data.par_order['r1'])*data.nvisit]
		self.r2 = params[data.par_order['r2']*data.nvisit:(1 + data.par_order['r2'])*data.nvisit]
		self.r3 = params[data.par_order['r3']*data.nvisit:(1 + data.par_order['r3'])*data.nvisit]
		self.scale = params[data.par_order['scale']*data.nvisit:(1 + data.par_order['scale'])*data.nvisit]
		self.amp1 = params[data.par_order['amp1']*data.nvisit:(1 + data.par_order['amp1'])*data.nvisit]
		self.theta1 = params[data.par_order['theta1']*data.nvisit:(1 + data.par_order['theta1'])*data.nvisit]
		self.amp2 = params[data.par_order['amp2']*data.nvisit:(1 + data.par_order['amp2'])*data.nvisit]
		self.theta2 = params[data.par_order['theta2']*data.nvisit:(1 + data.par_order['theta2'])*data.nvisit]
		self.E_s = params[data.par_order['E_s']*data.nvisit:(1 + data.par_order['E_s'])*data.nvisit]
		self.E_f = params[data.par_order['E_f']*data.nvisit:(1 + data.par_order['E_f'])*data.nvisit]
		self.eta_s = params[data.par_order['eta_s']*data.nvisit:(1 + data.par_order['eta_s'])*data.nvisit]
		self.eta_f = params[data.par_order['eta_f']*data.nvisit:(1 + data.par_order['eta_f'])*data.nvisit]
		self.tau_s = params[data.par_order['tau_s']*data.nvisit:(1 + data.par_order['tau_s'])*data.nvisit]
		self.tau_f = params[data.par_order['tau_f']*data.nvisit:(1 + data.par_order['tau_f'])*data.nvisit]
		self.E_s0 = params[data.par_order['E_s0']*data.nvisit:(1 + data.par_order['E_s0'])*data.nvisit]
		self.E_f0 = params[data.par_order['E_f0']*data.nvisit:(1 + data.par_order['E_f0'])*data.nvisit]
		self.f = params[data.par_order['f']*data.nvisit:(1 + data.par_order['f'])*data.nvisit]
		self.A1 = params[data.par_order['A1']*data.nvisit:(1 + data.par_order['A1'])*data.nvisit]
		self.P1 = params[data.par_order['P1']*data.nvisit:(1 + data.par_order['P1'])*data.nvisit]
		self.phi1 = params[data.par_order['phi1']*data.nvisit:(1 + data.par_order['phi1'])*data.nvisit]
		self.A2 = params[data.par_order['A2']*data.nvisit:(1 + data.par_order['A2'])*data.nvisit]
		self.P2 = params[data.par_order['P2']*data.nvisit:(1 + data.par_order['P2'])*data.nvisit]
		self.phi2 = params[data.par_order['phi2']*data.nvisit:(1 + data.par_order['phi2'])*data.nvisit]

def PrintParams(m, data): 
	print "per\t", m.params[data.par_order['per']*data.nvisit:(1 + data.par_order['per'])*data.nvisit]
	print "t0\t", m.params[data.par_order['t0']*data.nvisit:(1 + data.par_order['t0'])*data.nvisit]
	print "t_s\t", m.params[data.par_order['t_secondary']*data.nvisit:(1 + data.par_order['t_secondary'])*data.nvisit]
	print "w\t", m.params[data.par_order['w']*data.nvisit:(1 + data.par_order['w'])*data.nvisit]
	print "a\t", m.params[data.par_order['a']*data.nvisit:(1 + data.par_order['a'])*data.nvisit]
	print "inc\t", m.params[data.par_order['inc']*data.nvisit:(1 + data.par_order['inc'])*data.nvisit]
	print "rp\t", m.params[data.par_order['rp']*data.nvisit:(1 + data.par_order['rp'])*data.nvisit]
	print "fp\t", m.params[data.par_order['fp']*data.nvisit:(1 + data.par_order['fp'])*data.nvisit]
	print "u1\t", m.params[data.par_order['u1']*data.nvisit:(1 + data.par_order['u1'])*data.nvisit]
	print "u2\t", m.params[data.par_order['u2']*data.nvisit:(1 + data.par_order['u2'])*data.nvisit]
	print "ecc\t", m.params[data.par_order['ecc']*data.nvisit:(1 + data.par_order['ecc'])*data.nvisit]
	print "c\t", m.params[data.par_order['c']*data.nvisit:(1 + data.par_order['c'])*data.nvisit]
	print "v\t", m.params[data.par_order['v']*data.nvisit:(1 + data.par_order['v'])*data.nvisit]
	print "v2\t",  m.params[data.par_order['v2']*data.nvisit:(1 + data.par_order['v2'])*data.nvisit]
	print "r1\t", m.params[data.par_order['r1']*data.nvisit:(1 + data.par_order['r1'])*data.nvisit]
	print "r2\t", m.params[data.par_order['r2']*data.nvisit:(1 + data.par_order['r2'])*data.nvisit]
	print "r3\t", m.params[data.par_order['r3']*data.nvisit:(1 + data.par_order['r3'])*data.nvisit]
	print "scale\t", m.params[data.par_order['scale']*data.nvisit:(1 + data.par_order['scale'])*data.nvisit]
	print "amp1\t", m.params[data.par_order['amp1']*data.nvisit:(1 + data.par_order['amp1'])*data.nvisit]
	print "theta1\t", m.params[data.par_order['theta1']*data.nvisit:(1 + data.par_order['theta1'])*data.nvisit]
	print "amp2\t", m.params[data.par_order['amp2']*data.nvisit:(1 + data.par_order['amp2'])*data.nvisit]
	print "theta2\t", m.params[data.par_order['theta2']*data.nvisit:(1 + data.par_order['theta2'])*data.nvisit]
	print "E_s\t", m.params[data.par_order['E_s']*data.nvisit:(1 + data.par_order['E_s'])*data.nvisit]
	print "E_f\t", m.params[data.par_order['E_f']*data.nvisit:(1 + data.par_order['E_f'])*data.nvisit]
	print "eta_s\t", m.params[data.par_order['eta_s']*data.nvisit:(1 + data.par_order['eta_s'])*data.nvisit]
	print "eta_f\t", m.params[data.par_order['eta_f']*data.nvisit:(1 + data.par_order['eta_f'])*data.nvisit]
	print "tau_s\t", m.params[data.par_order['tau_s']*data.nvisit:(1 + data.par_order['tau_s'])*data.nvisit]
	print "tau_f\t", m.params[data.par_order['tau_f']*data.nvisit:(1 + data.par_order['tau_f'])*data.nvisit]
	print "E_s0\t", m.params[data.par_order['E_s0']*data.nvisit:(1 + data.par_order['E_s0'])*data.nvisit]
	print "E_f0\t", m.params[data.par_order['E_f0']*data.nvisit:(1 + data.par_order['E_f0'])*data.nvisit]
	print "f\t", m.params[data.par_order['f']*data.nvisit:(1 + data.par_order['f'])*data.nvisit]
	print "A1\t", m.params[data.par_order['A1']*data.nvisit:(1 + data.par_order['A1'])*data.nvisit]
	print "P1\t", m.params[data.par_order['P1']*data.nvisit:(1 + data.par_order['P1'])*data.nvisit]
	print "phi1\t", m.params[data.par_order['phi1']*data.nvisit:(1 + data.par_order['phi1'])*data.nvisit]
	print "A2\t", m.params[data.par_order['A2']*data.nvisit:(1 + data.par_order['A2'])*data.nvisit]
	print "P2\t", m.params[data.par_order['P2']*data.nvisit:(1 + data.par_order['P2'])*data.nvisit]
	print "phi2\t", m.params[data.par_order['phi2']*data.nvisit:(1 + data.par_order['phi2'])*data.nvisit]
	print "perror", m.perror

class Model:
	"""
	doc
	"""
	def __init__(self, params, data, flags):
		p = FormatParams(params, data)
		self.lc = np.zeros(len(data.time))				#full model light curve (with systematics)
		self.lc_nosys = np.zeros(len(data.time))			#full model light curve (no systematics)
		self.transit_model = np.zeros(len(data.time))			#transit model (relative flux; no systematics)
		self.eclipse_model = np.zeros(len(data.time))			#eclipse model (relative flux; no systematics)
		self.phase_model = np.zeros(len(data.time))			#phase curve model (relative flux; no systematics; includes eclipse)
		self.phase_model_no_eclipse = np.zeros(len(data.time))		#phase variation model (relative flux; no systematics; no eclipse))
		self.data_corr = np.zeros(len(data.time))			#data with the odd/even effect and orbit-long ramps removed
		self.phase = np.zeros(len(data.time))				#orbital phase	(defined in model because it depends on the period and ephemeris)
		self.sys = np.zeros(len(data.time))				#systematics showing visit-long trends only
		self.all_sys = np.zeros(len(data.time))				#all systematics 
		self.data_normalized = np.zeros(len(data.time))

		for i in range(data.nvisit):
			ind = data.vis_num == i
			self.phase[ind] = (data.time[ind] - p.t0[i])/p.per[i] - np.floor((data.time[ind] - p.t0[i])/p.per[i])
			if data.lc_type == "eclipse": 
				self.lc[ind] = get_elc(data.time[ind], p, data, i) 	
				self.lc_nosys[ind] = get_elc(data.time[ind], p, data, i) 	
			elif data.lc_type == "transit": 
				self.lc[ind] = get_tlc(data.time[ind], p, data, i) 
				self.lc_nosys[ind] = get_tlc(data.time[ind], p, data, i) 
			elif data.lc_type == "phase_curve":
				self.eclipse_model[ind] = get_elc(data.time[ind], p, data, i) 	
				self.transit_model[ind] = get_tlc(data.time[ind], p, data, i) 
				self.phase_model[ind] = get_phaselc(data.time[ind], p, data, i)
				self.lc[ind] = self.transit_model[ind] + (self.eclipse_model[ind] - 1.)*self.phase_model[ind]
				self.lc_nosys[ind] = self.transit_model[ind] + (self.eclipse_model[ind] - 1.)*self.phase_model[ind]
			else: assert False, "Unknown option; supported light curve types are 'transit', 'eclipse', and 'phase_curve'"
			self.phase_hr = np.linspace(self.phase.min()-0.05, self.phase.max()+0.05, 1000)
			self.t_hr = self.phase_hr*p.per[0]+p.t0[0]
			self.lc_hr = get_tlc(self.t_hr, p, data, 0)

			if flags['divide-white'] == False:
				S = np.ones_like(self.transit_model[ind])+p.scale[i]*data.scan_direction[ind]
				D = np.ones_like(self.transit_model[ind])+p.r3[i]*data.t_delay[ind]
				self.all_sys[ind] = data.flux[ind]/self.lc[ind]
				self.lc[ind] *= (p.c[i]*S + p.v[i]*data.t_vis[ind] + p.v2[i]*data.t_vis[ind]**2)*(1.0-np.exp((-p.r1[i]*data.t_orb[ind]-p.r2[i]-D)))*(p.A1[i]*np.sin(2.*np.pi*(data.t_vis[ind] - p.phi1[i])/p.P1[i]))
				self.data_corr[ind] = data.flux[ind]/((p.c[i]*S + p.v[i]*data.t_vis[ind] + p.v2[i]*data.t_vis[ind]**2)*(1.0-np.exp(-p.r1[i]*data.t_orb[ind]-p.r2[i]-D))*(p.A1[i]*np.sin(2.*np.pi*(data.t_vis[ind] - p.phi1[i])/p.P1[i])))
				self.sys[ind] = p.c[i]+ p.v[i]*data.t_vis[ind] + p.v2[i]*data.t_vis[ind]**2
				self.data_normalized[ind] = data.flux[ind]/(1.-np.exp((-p.r1[i]*data.t_orb[ind]-p.r2[i]-D))) - self.transit_model[ind]*p.c[i]*(S-1.)
			#	E_s = E(data.t_vis[ind], p.eta_s[i], p.f[i], p.E_s[i], p.tau_s[i], p.E_s0[i])
			#	E_f = E(data.t_vis[ind], p.eta_f[i], p.f[i], p.E_f[i], p.tau_f[i], p.E_f0[i])
			#	self.lc[ind] *= p.c[i]*(p.f[i] - (p.eta_s[i]*p.f[i] - E_s*(p.eta_s[i]*p.f[i]/p.E_s[i] + 1./p.tau_s[i])*np.exp(-data.t_vis[ind]*(p.eta_s[i]*p.f[i]/p.E_s[i] + 1./p.tau_s[i]))) - E_f*(p.eta_f[i]*p.f[i]/p.E_f[i] + 1./p.tau_f[i])*np.exp(-data.t_vis[ind]*(p.eta_f[i]*p.f[i]/p.E_f[i] + 1./p.tau_f[i])))
			#	self.data_corr[ind] = data.flux[ind]/(p.c[i]*(p.f[i] - (p.eta_s[i]*p.f[i] - E_s*(p.eta_s[i]*p.f[i]/p.E_s[i] + 1./p.tau_s[i])*np.exp(-data.t_vis[ind]*(p.eta_s[i]*p.f[i]/p.E_s[i] + 1./p.tau_s[i]))) - E_f*(p.eta_f[i]*p.f[i]/p.E_f[i] + 1./p.tau_f[i])*np.exp(-data.t_vis[ind]*(p.eta_f[i]*p.f[i]/p.E_f[i] + 1./p.tau_f[i]))))
				#self.sys[ind] = p.c[i]+ p.v[i]*data.t_vis[ind] + p.v2[i]*data.t_vis[ind]**2
				#self.data_normalized[ind] = data.flux[ind]/(1.-np.exp((-p.r1[i]*data.t_orb[ind]-p.r2[i]-D))) - self.transit_model[ind]*p.c[i]*(S-1.)
			else:
				self.lc[ind] *= (p.c[i] + p.v[i]*data.t_vis[ind] + p.v2[i]*data.t_vis[ind]**2)*data.all_sys[ind]
				self.data_corr[ind] = data.flux[ind]/((p.c[i] + p.v[i]*data.t_vis[ind] + p.v2[i]*data.t_vis[ind]**2)*data.all_sys[ind])

		#plt.plot(data.time, self.lc, 'xk')
		#plt.show()

		self.resid = data.flux - self.lc 
		self.chi2 = np.sum((self.resid/data.err)**2)		
		self.chi2red = self.chi2/data.dof
	#	ind = data.orb_num!=2
	#	print "calculating rms without orbit 2"
	#	self.rms = 1.0e6*np.sqrt(np.mean((self.resid[ind]/data.flux[ind])**2))
	#	self.rms_predicted = 1.0e6*np.sqrt(np.mean((data.err[ind]/data.flux[ind])**2))
		self.rms = 1.0e6*np.sqrt(np.mean((self.resid/data.flux)**2))
		self.rms_predicted = 1.0e6*np.sqrt(np.mean((data.err/data.flux)**2))
		self.data = data
		self.ln_likelihood = -0.5*(np.sum((self.resid/data.err)**2 + np.log(2.0*np.pi*(data.err)**2)))
		self.bic = -2.*self.ln_likelihood + data.nfree_param*np.log(data.npoints)


def plot(params, data, flags, obs_par, plot_sys=False):
	p = FormatParams(params, data)
	m = Model(params,data, flags)
	sns.set_palette("muted")
	palette = sns.color_palette("muted", m.data.nvisit)
	if plot_sys==False:
		for i in range(m.data.nvisit):
			ind = m.data.vis_num == i
			plt.subplot(211)
			if data.lc_type == "transit":
				temp_phase = np.copy(m.phase[ind])
				temp_phase[temp_phase>0.5] -= 1.0
			 	plt.plot(temp_phase, m.data_corr[ind], marker='o', linestyle="None", markersize=5)
			else: plt.plot(m.phase[ind], m.data_corr[ind], marker='o', linestyle="None", markersize=5)
			plt.subplot(212)
			if data.lc_type == "transit": plt.plot(temp_phase, 1.0e6*m.resid[ind]/m.data.flux[ind], marker='o', linestyle="None", markersize=5)
			else: plt.plot(m.phase[ind], 1.0e6*m.resid[ind]/m.data.flux[ind], marker='o', linestyle="None", markersize=5)
		plt.subplot(211)
		phase_hr = np.linspace(m.phase.min()-0.05, m.phase.max()+0.05, 1000)
		t_hr = phase_hr*p.per[0]+p.t0[0]

		if data.lc_type == "transit":
			transit_model_hr = get_tlc(t_hr, p, data, 0) 
			lc_hr = transit_model_hr
		elif data.lc_type == "eclipse":
			eclipse_model_hr = get_elc(t_hr, p, data, 0) 
			lc_hr = eclipse_model_hr
		elif data.lc_type == "phase_curve":
			transit_model_hr = get_tlc(t_hr, p, data, 0) 
			eclipse_model_hr = get_elc(t_hr, p, data, 0) 
			phase_model_hr = get_phaselc(t_hr, p, data, 0)
			lc_hr = transit_model_hr + (eclipse_model_hr-1.)*phase_model_hr
		else: assert False, "Unknown option; supported light curve types are 'transit', 'eclipse', and 'phase_curve'"
		
		if data.lc_type == "transit": 
			phase_hr[phase_hr > 0.5] -= 1.0
			ind = np.argsort(phase_hr)
			phase_hr = phase_hr[ind]
			lc_hr = lc_hr[ind]
			
		plt.plot(phase_hr, lc_hr, color='0.2', zorder=1)
		plt.ylabel("Relative flux")
		ax = plt.gca()
		ax.text(0,1,'obs, exp rms (ppm); chi:\n '+'{0:d}'.format(int(m.rms))+", "+'{0:d}'.format(int(m.rms_predicted))+", "+'{0:0.2f}'.format(m.chi2red),verticalalignment='top',horizontalalignment='left',transform=ax.transAxes, fontsize=10) #LK comment!!
		if flags['fit-white']: plt.title("Fit to white light curve")
		else: plt.title("Fit to {0:0.2f}".format(m.data.wavelength)+" micron channel")
		delta = 3.0*m.rms/1.0e6
		plt.ylim((lc_hr.min() - delta , lc_hr.max()+delta))
		if data.lc_type == "phase_curve": plt.ylim((1.0-delta, lc_hr.max()+delta))
		plt.xlim((phase_hr.min(), phase_hr.max()))
		if obs_par['lc_type'] == 'transit': plt.xlim(-0.03, 0.03)
		ax = plt.gca()
		ax.ticklabel_format(useOffset=False)

		plt.subplot(212)
		plt.axhline(0, zorder=1, color='0.2', linestyle='dashed')
		plt.ylabel("Residuals (ppm)")
		plt.xlabel("Orbital phase")
		plt.xlim((phase_hr.min(), phase_hr.max()))
		if obs_par['lc_type'] == 'transit': plt.xlim(-0.03, 0.03)
		plt.show()	#LK comment
		#plt.savefig("temp.pdf")

		"""plt.hist(m.resid/m.data.flux*1e6, normed=True)

		xx = np.linspace(-2000, 2000, 100)
		plt.plot(xx, norm.pdf(xx, 0, int(m.rms_predicted)))

		plt.show()"""
		#make rms vs. bin size plot

	elif plot_sys == True:
		#offset = np.mean(m.sys)*0.001
		print "FIXME: this option isn't implemented for phase curves"
		offset = 0.
		for i in range(m.data.nvisit):
			ind = m.data.vis_num == i
			plt.subplot(211)
			plt.plot(m.phase[ind], m.data_normalized[ind]+offset*i, marker='o', linestyle="None", markersize=5, color=palette[i])
			plt.plot(m.phase[ind], m.sys[ind]*m.transit_model[ind]+offset*i, color=palette[i], zorder=1)
			plt.plot(m.phase[ind], m.sys[ind]+offset*i, color=palette[i], zorder=1, linestyle='dashed')
			plt.ylabel("Corrected flux")
			plt.title("Flux corrected for orbit-long trends and scan direction")
			plt.subplot(212)
			plt.plot(m.phase[ind], m.resid[ind], marker='o', linestyle="None", markersize=5, color=palette[i])
			plt.xlabel("Orbital phase")
			plt.ylabel("Residuals (e-)")
		plt.subplot(211)
		plt.plot(m.phase[0], m.sys[0]*m.transit_model[0], color='0.5', zorder=1, label="with transit")
		plt.plot(m.phase[0], m.sys[0], color='0.5', zorder=1, linestyle='dashed', label="transit removed")
		plt.legend()
		plt.show()
			
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

def plot_data(data):
	palette = sns.color_palette("husl", data.nvisit)
	for i in range(data.nvisit): 	
		ind = data.vis_num==i
		plt.subplot((data.nvisit)*100+10+i+1)
		plt.plot(data.t_vis[ind]*24., data.flux[ind], marker='o', markersize=4.5, linestyle="none", label = "Visit {0}".format(i), color= palette[i])
		plt.xlim(((data.t_vis.min()-0.02)*24., (data.t_vis.max()+0.05)*24.))
		plt.ylim((0.998*data.flux.min(), 1.002*data.flux.max()))
		plt.legend()
	plt.xlabel("Time after visit start (hours)")
	plt.ylabel("Flux (e-)")
	plt.tight_layout()
	plt.show()	


def weighted_mean(data, err):            #calculates the weighted mean for data points data with std devs. err
	ind = err != 0.0
        weights = 1.0/err[ind]**2
        mu = np.sum(data[ind]*weights)/np.sum(weights)
        var = 1.0/np.sum(weights)
        return [mu, np.sqrt(var)]                

def get_tlc(t, p, data, v_num):
	bat_params = batman.TransitParams()
	bat_params.t0 = p.t0[v_num]
	bat_params.t_secondary = p.t_secondary[v_num]
	bat_params.per = p.per[v_num]
	bat_params.rp = p.rp[v_num]
	bat_params.a = p.a[v_num]
	bat_params.inc = p.inc[v_num]
	bat_params.ecc = p.ecc[v_num]
	bat_params.w = p.w[v_num]
	bat_params.u = [p.u1[v_num], p.u2[v_num]]
	bat_params.limb_dark = "quadratic"		#FIXME - specify this value in one of the config files
	
	m = batman.TransitModel(bat_params, t, supersample_factor=3, exp_time = data.exp_time/24./60./60.)
	return m.light_curve(bat_params)

def get_elc(t, p, data, v_num):
	bat_params = batman.TransitParams()
	bat_params.t0 = p.t0[v_num]
	bat_params.t_secondary = p.t_secondary[v_num]
	bat_params.per = p.per[v_num]
	bat_params.rp = p.rp[v_num]
	bat_params.a = p.a[v_num]
	bat_params.inc = p.inc[v_num]
	bat_params.ecc = p.ecc[v_num]
	bat_params.w = p.w[v_num]
	bat_params.u = [p.u1[v_num], p.u2[v_num]]
	bat_params.fp = p.fp[v_num] 
	bat_params.limb_dark = "quadratic"		#FIXME - specify this value in one of the config files
	
	m = batman.TransitModel(bat_params, t, supersample_factor=3, exp_time = data.exp_time/24./60./60., transittype="secondary")
	return m.light_curve(bat_params)

def get_phaselc(t, p, data, v_num):
	return 1.+p.amp1[v_num]*np.cos(2.*np.pi*(t-p.theta1[v_num])/p.per[v_num]) + p.amp2[v_num]*np.cos(4.*np.pi*(t-p.theta2[v_num])/p.per[v_num])

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
	
	if flags['plot-raw-data']: plot_data(data)
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

	if flags['show-plot']: plot(m.params, data, flags, obs_par, plot_sys=flags['plot-sys'])
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
