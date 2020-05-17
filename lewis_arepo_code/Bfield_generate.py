import numpy as np
import scipy.special as sc
import  scipy.optimize as sp 
import scipy.interpolate as interp

def spectrum(theta,kstar,Rm_crit,boxsize,N,rho,v_turb):
	'''turbulent dynamo saturated B field energy spectrum from Schober(2015)
	set k (2pi/lambda) limits as lambda=boxsize to lambda=boxsize/(N-1), where N is the 
	number of pixels you want in the real space image'''
	kmin=2*np.pi/boxsize
	kmax=2*np.pi/(boxsize/(N-1))
	kL=kmin
	
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


def discrete_modes(k,M,boxsize_cgs):
	'''discrete fourier transfer requires integral over k where
	k is the number of cycles over the box length in pixels, from 0 to N-1 where N 
	is the number of pixels. The k from spectrum() is a wavenumber 2pi/lambda (cgs). 
	This function converts wavenumbers into new k.'''
	k_new=k*boxsize_cgs/(2*np.pi)

	#k_new inversely proportional to k_old, needs rearranging 
	#M_new=M[np.argsort(k_new)]
	#k_new=k_new[np.argsort(k_new)]
	return k_new


def modes(k,M,N):
	'''generate 3D wavevectors, aplitudes and phases, removes divergent modes
	-based on Lomax(2015)'''
	
	#first fit a spline to power spectrum
	M_function=interp.interp1d(k,M,bounds_error=False,fill_value=0)
	
	#3D wavevectors have a 3D amplitude for each k=(kx,ky,kz)	

	N=int(N)
	#need space for B field distribution
	Bx=np.zeros((N,N,N),dtype=np.complex128)
	By=np.zeros((N,N,N),dtype=np.complex128)
	Bz=np.zeros((N,N,N),dtype=np.complex128)

        #values for navigating k space
	ks=np.linspace(0,N-1,N)

	#for each coordinate in 3D k space, assign a 3D amplitude and phase 
	print('creating k space')
	for X in range(N):
		for Y in range(N):
			for Z in range(N):
				#figure out where we are on the power spectrum
				K=np.sqrt( (ks[X])**2 + (ks[Y])**2 + (ks[Z])**2)
				
				P=M_function(K)
				#sample amplitude from Gaussian
				Ax=np.sqrt(P)*np.random.normal(0,1)
				Ay=np.sqrt(P)*np.random.normal(0,1)
				Az=np.sqrt(P)*np.random.normal(0,1)
				
				#phase from random distribution
				phix=2*np.pi*np.random.uniform(0, 1)
				phiy=2*np.pi*np.random.uniform(0, 1)
				phiz=2*np.pi*np.random.uniform(0, 1)
	
				#remove compressive modes (no div B)
				#dot product between k and A
				dotx=(Ax*(ks[X]))
				doty=(Ay*(ks[Y]))
				dotz=(Az*(ks[Z]))
				dot_product=dotx+doty+dotz
				#refer to eq 2.6/2.7 Lomax(2015)
				Ax=Ax-(ks[X])*dot_product
				Ay=Ay-(ks[Y])*dot_product
				Az=Az-(ks[Z])*dot_product

				#use phase to split into real/imaginary
				Bx[X,Y,Z]=Ax*complex(np.cos(phix),np.sin(phix))
				By[X,Y,Z]=Ay*complex(np.cos(phiy),np.sin(phiy))
				Bz[X,Y,Z]=Az*complex(np.cos(phiz),np.sin(phiz))
	
	print('performing transform')	
	#reverse FFT into real space
	bx=np.fft.ifftn(Bx,np.array([N,N,N]))
	by=np.fft.ifftn(By,np.array([N,N,N]))
	bz=np.fft.ifftn(Bz,np.array([N,N,N]))


	return Bx,By,Bz,bx,by,bz
