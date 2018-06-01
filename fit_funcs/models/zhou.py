import numpy as np

def E(t, eta, f, Etot, tau, E0):
	return eta*f/(eta*f/Etot + 1./tau) + (E0 - eta*f/(eta*f/Etot + 1./tau))*np.exp(-(eta*f/Etot +1./tau)*t)

