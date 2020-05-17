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
	kz=np.linspace(0,N-1,N)
	ky,kx=np.meshgrid(kz,kz)

	print('creating k space')
	#for each coordinate in 3D k space, assign a 3D amplitude and phase
	#split k space sube into xy sheet stack
	for z in range(N): 
	
		if z==int(N/2):
			print('50%')
		
		#figure out where on the power spectrum we are 
		K=np.sqrt( kx**2 + ky**2 + (kz[z])**2)
		P=M_function(K)

		#assign 3D amplitude from guassian distribution
		Ax=np.sqrt(P)*np.random.normal(0,1,(N,N))
		Ay=np.sqrt(P)*np.random.normal(0,1,(N,N))
		Az=np.sqrt(P)*np.random.normal(0,1,(N,N))

		#assign 3D phase from uniform distribution
		phix=2*np.pi*np.random.uniform(0, 1,(N,N))
		phiy=2*np.pi*np.random.uniform(0, 1,(N,N))
		phiz=2*np.pi*np.random.uniform(0, 1,(N,N))	

		#remove the divergent modes 
		#refer to eq 2.6/2.7 from Lomax(2015)
		dotx=(Ax*(kx))
		doty=(Ay*(ky))
		dotz=(Az*(kz[z]))
		dot_product=dotx+doty+dotz
		Ax=Ax-(kx)*dot_product
		Ay=Ay-(ky)*dot_product
		Az=Az-(kz[z])*dot_product

		#split the signal into real/imaginary parts based on phase 
		Bx[:,:,z]=Ax*(np.cos(phix)+np.sin(phix)*1j)
		By[:,:,z]=Ay*(np.cos(phiy)+np.sin(phiy)*1j)
		Bz[:,:,z]=Az*(np.cos(phiz)+np.sin(phiz)*1j)
	
	print('performing transform')	
	#reverse FFT into real space
	bx=np.fft.ifftn(Bx,np.array([N,N,N]))
	by=np.fft.ifftn(By,np.array([N,N,N]))
	bz=np.fft.ifftn(Bz,np.array([N,N,N]))


	return Bx,By,Bz,bx,by,bz
