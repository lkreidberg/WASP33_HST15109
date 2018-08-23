import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt

h = 6.626e-34   #J/s
c = 3.0e8       #m/s
kb = 1.4e-23     #J/K

def get_significance(chisq, dof):
	alpha = (1. - st.chi2.cdf(chisq, dof))/2.
	z = st.norm.ppf(1.-alpha)
	return z


def blackbody(l,T): 
	return 2*h*c**2/(l**5*(np.exp(h*c/(l*kb*T))-1))


def best_fit_bb(w, y, e, Tstar, rprs):
	Ts = np.linspace(1800, 4000, 100)
	chibest = 10000.
	Tbest = 0.	
	w = np.array(w)
	for T in Ts:
		model = 1.0e6*blackbody(w*1.0e-6, T)/blackbody(w*1.0e-6, Tstar)*rprs**2
		chi2 = np.sum((y - model)**2/e**2)
		if chi2 < chibest: 
			chibest, Tbest = chi2, T
	#		print("chi2, Tbest", chibest/(len(w)-1.)), Tbest
	waves_hires = np.linspace(0.7, 1.2, 1000)
	#return waves_hires, 1.0e6*blackbody(waves_hires*1.0e-6, Tbest)/blackbody(waves_hires*1.0e-6, Tstar)*rprs**2, Tbest, get_significance(chibest,(len(w)-1))
	return waves_hires, 1.0e6*blackbody(waves_hires*1.0e-6, Tbest)/blackbody(waves_hires*1.0e-6, Tstar)*rprs**2, Tbest, chibest/(len(w)-1)


#files = ["../analysis/fit_2018_08_23_11:16.txt"]
#files = ["../analysis/fit_2018_08_23_11:16_cut.txt"]
#files = ["../analysis/fit_2018_08_23_16:36.txt"]
files = ["../analysis/fit_2018_08_23_16:55.txt"]        #analytic model
Tstars = [7400.] 
rps = [0.103] 

for i, f in enumerate(files):
	plt.figure(figsize = (5,2))
	Tstar = Tstars[i]
	rprs = rps[i]
	d = np.genfromtxt(f)
        #w, y, e = d[:,0], d[:,1]*1e6, d[:,2]*np.sqrt(d[:,5])*1e6
        w, y, e = d[:,0], d[:,1]*1e6, d[:,2]*1e6

	xm, ym, Tbest, chi2  = best_fit_bb(w, y, e, Tstar, rprs)

	plt.plot(xm, ym, label = str(Tbest)+ " blackbody", color='0.7')
	plt.errorbar(w, y, e, fmt = 'ok')
	#plt.legend()

	print f, Tbest, chi2
	#plt.savefig("bb_fit.pdf")
	plt.show()
