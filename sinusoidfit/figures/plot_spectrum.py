import numpy as np
import matplotlib.pyplot as plt
from astropy.convolution import convolve, Gaussian1DKernel
import seaborn as sns
import scipy.stats as st

sns.set_context("talk", font_scale=1.5)
sns.set_style("white")
sns.set_style("ticks", {"xtick.direction":"in", "ytick.direction":"in"})

h = 6.626e-34   #J/s
c = 3.0e8       #m/s
kb = 1.4e-23     #J/K


def boxcar_smooth(x, nsmooth):
	n = len(x)
	for i in np.arange(0, n):
		lower = np.max(np.array([0, i - int(nsmooth/2)]))
		upper = np.min(np.array([n-1, i + int(nsmooth/2)]))
		x[i] = 	np.mean(x[lower:upper])
	return x


def get_significance(chisq, dof):
	alpha = (1. - st.chi2.cdf(chisq, dof))/2.
	z = st.norm.ppf(1.-alpha)
	return z


def blackbody(l,T): 
	return 2*h*c**2/(l**5*(np.exp(h*c/(l*kb*T))-1))


def best_fit_bb(w, y, e, Tstar, rprs):
	Ts = np.linspace(2800, 3200, 1000)
	chibest = 10000.
	Tbest = 0.	
	w = np.array(w)
	for T in Ts:
		model = blackbody(w*1.0e-6, T)/blackbody(w*1.0e-6, Tstar)*rprs**2
		chi2 = np.sum((y - model)**2/e**2)
		if chi2 < chibest: 
			chibest, Tbest = chi2, T
	#		print("chi2, Tbest", chibest/(len(w)-1.)), Tbest
	waves_hires = np.linspace(0.5, 1.7, 100)
	return waves_hires, blackbody(waves_hires*1.0e-6, Tbest)/blackbody(waves_hires*1.0e-6, Tstar)*rprs**2, Tbest, get_significance(chibest,(len(w)-1))



g = Gaussian1DKernel(stddev=1.1)

files  = ["wasp33b_from_Madhu/spec.dat", "wasp33b_from_Madhu/spec_notio.dat"]
labels  = ["best fit model", "best fit with no TiO"]
colors = ['blue', 'red']
linestyles = ['dotted', 'dashed']

plt.figure(figsize = (12,6))


"""for i, f in enumerate(files):
	d = np.genfromtxt(f, skip_header = 2)
	plt.plot(d[:,1]*1.e4, convolve(d[:,4], g, boundary = 'extend'), label = labels[i], color = '0.5', alpha = 0.5, linestyle = linestyles[i])"""

s = np.genfromtxt("w33_data_haynes.txt")
plt.errorbar(s[:,0], s[:,1]/100., s[:,2]/100., fmt = 'ow', markersize = 4, ecolor = 'k', markeredgecolor = 'k', markeredgewidth = 1., linewidth = 1., linestyle='none', zorder=100, label="G141 data (Haynes)")

#s = np.genfromtxt("../../WASP33_HST12495/vis_combined/mcmc_11bins_mr_fixsine/emission_spectrum.txt")
s = np.genfromtxt("../../WASP33_HST12495/vis_combined/fit_2018_09_12_13:58.txt")
#offset = 0.0001
offset = 0.0
plt.errorbar(s[:,0], s[:,1] - offset, s[:,2], marker='.', color='r', linestyle='none', zorder=100, label="G141 data (Kreidberg)")


d = np.genfromtxt("../mcmc_output_ackbar_fixamps/emission_spectrum.txt")
plt.errorbar(d[:,0], d[:,1], d[:,2], fmt = '.b', zorder=100, label = "G102 data")

xm, ym, Tbest, chi2  = best_fit_bb(d[:,0], d[:,1], d[:,2], 7400, 0.1106)
print "best fit T", Tbest
plt.plot(xm, ym, color='0.5',  label = 'blackbody fit to G102', alpha = 0.5, linestyle='dotted', zorder = 0.5)


plt.ylim(0.0, 1.8e-3)
plt.xlim(0.75, 1.7)


#plt.gca().annotate('TiO features', xy=(1.05, 0.0008), xytext=(0.8, 0.0005), arrowprops=dict(facecolor='black', shrink=0.05),)
#plt.gca().annotate('', xy=(0.97, 0.00095), xytext=(0.95, 0.0006), arrowprops=dict(facecolor='black', shrink=0.05),)


plt.tight_layout()
plt.xlabel("Wavelength (microns)")
plt.ylabel("Planet-to-star flux")
plt.legend(loc = "lower right", frameon=True, fontsize=16)

plt.savefig("w33_models.pdf")
plt.show() 
