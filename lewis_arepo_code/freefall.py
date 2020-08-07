import numpy as np
import astropy.constants as ap

print('radius: 5.770217077343649e+18')
print('Mass: 5.415932109255394e+36')


def t(M,R):
	V=4/3*np.pi*R**3
	rho=M/V
	t_ff=np.sqrt(3*np.pi/(32*rho*ap.G.cgs.value))
	return t_ff
