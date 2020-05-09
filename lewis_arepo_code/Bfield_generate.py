import numpy as np
import scipy.special as sc

def spectrum(theta,kmin,kmax,kstar,kL,Rm_crit,rho,v_turb):
	gamma = (304*theta + 163) / 60
	logk=np.linspace(np.log10(kmin),np.log10(kmax),1000)
	k=10**logk
	dk=np.zeros_like(k)
	for i in range(1000-1):
		dk[i]=k[i+1]-k[i]
	M=k**(3/2)*sc.kn(1,k/kstar)
	alpha= 3/4 *gamma *(kL/kstar)**(2*theta)*rho*v_turb**2  /sum(M*dk)
	M_norm=alpha*M
	return M_norm,k
	

