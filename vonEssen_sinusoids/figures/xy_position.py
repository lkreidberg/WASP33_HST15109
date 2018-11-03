import numpy as np
import matplotlib.pyplot as plt

d = np.genfromtxt("../extracted_lc/08_25_12_53/diagnostics.txt")


plt.subplot(211)
plt.title("X position")
plt.plot(d[:,1], d[:,4], '.k')

plt.subplot(212)
plt.title("Y position")
ind = d[:,5] == -10.
plt.plot(d[~ind,1], d[~ind,5], '.k')
plt.show()
