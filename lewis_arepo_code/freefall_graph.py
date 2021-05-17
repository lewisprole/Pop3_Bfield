import numpy as np
import matplotlib.pyplot as plt 
from scipy.stats import binned_statistic
import code_units 
import astropy.constants as ap

'''note to self:
wanted 2000 years of free-fall time
gives r~0.25 code units 
double to make safe r~0.5 code units 
therefore boxsize = 1 CU
'''

def rs_rho(a):
	mid=np.where(a.mass/a.rho==min(a.mass/a.rho))
	r=np.sqrt((a.x-a.x[mid])**2+(a.y-a.y[mid])**2+(a.z-a.z[mid])**2)
	r=r*code_units.d_cu
	bins=10**(np.linspace(np.log10(np.sort(r)[1]),np.log10(r.max()),100))
	rho,rs,z=binned_statistic(r,a.rho*code_units.rho_cu,bins=bins)
	return rho,rs[:-1]

def t_ff(rho):
	return np.sqrt(3*np.pi/(32*ap.G.cgs.value*rho))

def plotit(a):
	rho,rs=rs_rho(a)
	tff=t_ff(rho)
	plt.figure()
	plt.loglog(rs/code_units.d_cu,tff/(60*60*24*365)) #r in code units
