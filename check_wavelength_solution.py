import sys
sys.path.insert(0, './util')
import matplotlib.pyplot as plt
import numpy as np
import pyfits
import geometry102 as geo

d = pyfits.open("/Users/lkreidberg/Desktop/Data/WASP107_HST14916/idex01n5q_ima.fits")

plt.plot(d[1].data.sum(axis = 0))
plt.show()

refpix = np.genfromtxt("config/xrefyref.txt")

xref = refpix[0,1]
yref = refpix[0,2]
LTV1 = d[1].header['LTV1']
LTV2 = d[1].header['LTV2']
BEAMA_i = 41
BEAMA_f = 248

trace_i = int(yref + BEAMA_i + LTV1)
trace_f = yref + BEAMA_f + LTV1
print("trace_i, trace_f", trace_i, trace_f)

flux = (d[1].data[:,int(trace_i):266]).sum(axis = 0)

delx = 0.5 + np.arange(266 - trace_i) + BEAMA_i
disp = geo.dispersion(xref, yref)				#determines dispersion coefficients
print("disp", disp)
w = disp[0] + delx*disp[1]

print("w", w)
plt.plot(w/1.e4, flux/np.max(flux))

s = np.load("G102_sensitivity.npy")
plt.plot(s[:,0], s[:,1], color='r')

plt.show()

