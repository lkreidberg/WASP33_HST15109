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
	Ts = np.linspace(2800, 4000, 100)
	chibest = 10000.
	Tbest = 0.	
	w = np.array(w)
	for T in Ts:
		model = blackbody(w*1.0e-6, T)/blackbody(w*1.0e-6, Tstar)*rprs**2
		print np.median(model)
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

plt.figure(figsize = (12,6))


for i, f in enumerate(files):
	d = np.genfromtxt(f, skip_header = 2)

	plt.plot(d[:,1]*1.e4, convolve(d[:,4], g, boundary = 'extend'), label = labels[i], color = colors[i])
	#plt.plot(d[:,1]*1.e4, d[:,4], label = labels[i]) 
#	plt.plot(d[:,1]*1.e4, boxcar_smooth(d[:,4], 5), label = labels[i], color = colors[i])

s = np.genfromtxt("w33_data_haynes.txt")
offset = 0.00017
plt.errorbar(s[:,0], s[:,1]/100. - offset, s[:,2]/100., marker='o', color='0.5', linestyle='none', zorder=100, label="G141 data (Haynes et al. 2015)")

xm, ym, Tbest, chi2  = best_fit_bb(s[:,0], s[:,1]/100. - offset, s[:,2]/100., 7400, 0.1106)
print Tbest
plt.plot(xm, ym, color='0.4', linestyle='dotted', label = 'blackbody fit')


#d = np.genfromtxt("../analysis/fit_2018_08_23_10:45.txt")
d = np.genfromtxt("../analysis/fit_2018_08_23_11:16.txt")
plt.errorbar(d[:,0], d[:,1], d[:,2], fmt = 'ok', zorder=100, label = "G102 data")


#plt.ylim(0, 1.3e-3)
plt.ylim(0.2e-3, 1.5e-3)
plt.xlim(0.7, 1.7)


plt.gca().annotate('TiO features', xy=(1.05, 0.0008), xytext=(0.8, 0.0005), arrowprops=dict(facecolor='black', shrink=0.05),)
plt.gca().annotate('', xy=(0.97, 0.00095), xytext=(0.95, 0.0006), arrowprops=dict(facecolor='black', shrink=0.05),)


plt.tight_layout()
plt.xlabel("Wavelength (microns)")
plt.ylabel("Planet-to-star flux")
plt.legend(loc = "lower right", frameon=True, fontsize=16)

plt.savefig("w33_models.pdf")
plt.show() 
