import numpy as np
import astropy.constants as ap



G=ap.G.cgs.value
kb=ap.k_B.cgs.value
mp=ap.m_p.cgs.value


def B_strength(M,R,gamma0):
	'''calculates field strength for a gamma=(B energy/ Grav energy)'''
	U=3/5*G*M**2/R
	B_en=1/(8*np.pi)* 4/3 * np.pi * R**3 #magnetic energy/B^2
	B=np.sqrt(gamma0*U/B_en)
	return B

def jeans_length(T,rho,mu):
	'''calculates jeans length for temperature, density and mean weight'''
	return np.sqrt(kb*T/(G*rho*mp*mu))



