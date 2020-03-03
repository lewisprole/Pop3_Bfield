import numpy as np
import astropy.constants as ap

def B_strength(M,R,gamma0):
	G=ap.G.cgs.value
	U=3/5*G*M**2/R
	B_en=1/(8*np.pi)* 4/3 * np.pi * R**3 #magnetic energy/B^2
	B=np.sqrt(gamma0*U/B_en)
	return B

