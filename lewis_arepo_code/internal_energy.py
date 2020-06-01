import numpy as np


def int_en(N,T,mu):
	u=np.ones(N)
	R = 8.314e7
	gg = 6.672041e-8
	umass = 1.991e33
	udist = 1.000e17
	udens = umass/(udist)**3
	utime = np.sqrt((udist)**3/(gg*umass))
	uergg=(udist)**2/(utime)**2
    
	
	U = u*3.0/2.0*T*R/mu/uergg

	return U
#u=int_en(10000,10)


