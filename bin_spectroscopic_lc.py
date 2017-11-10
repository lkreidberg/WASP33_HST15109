#This code reads in the optimally extracted lightcurve and bins it into channels 5 pixels wide, following Berta '12
from numpy import *
from pylab import *
from astropy.io import ascii

def make_dict(table): return {x['parameter']: x['value'] for x in table}

def weighted_mean(data, err):				#calculates the weighted mean for data points data with variances err
	weights = 1.0/err
	mu = np.sum(data*weights)/np.sum(weights)
	var = 1.0/np.sum(weights)
	return [mu, var]				#returns weighted mean and variance

print "WARNING: THIS CODE DOES NOT INTERPOLATE AT THE EDGES OF THE BIN"

#reads in spectra
d = np.genfromtxt("extracted_lc/11_06_15_08/lc_spec.txt")

obs_par = make_dict(ascii.read("config/obs_par.txt", Reader=ascii.CommentedHeader))
nexp = int(obs_par['nexp'])			#number of exposures
npix = 207 				#width of spectrum in pixels (BEAMA_f - BEAMA_i) (sometimes, a little shorter because spectrum goes off edge of detector)
d = d.reshape(nexp, npix, -1)			#reshapes array by exposure

wave_bins = np.linspace(0.84, 1.1, 11)*1e4

#stores the indices corresponding to the wavelength range in each bin
wave_inds = []
for i in range(len(wave_bins)- 1): wave_inds.append((d[0,:,4] >= wave_bins[i])&(d[0,:,4] <= wave_bins[i+1]))

for i in range(len(wave_bins) - 1):
	wave = (wave_bins[i] + wave_bins[i+1])/2./1.e4
	outname = "speclc" + "{0:.3f}".format(wave)+".txt"
	outfile = open(outname, 'w')
	for j in range(nexp):
		time, phase, visnum, orbnum, scan = d[j,0,0], d[j,0,1], d[j,0,5], d[j,0,6], d[j,0,7]

		fluxes = d[j, wave_inds[i], 2]
		errs = d[j, wave_inds[i], 3]

		flux, variance = weighted_mean(fluxes, errs)		

		print>>outfile, phase, flux, variance, wave, 0., time, visnum, orbnum, scan


