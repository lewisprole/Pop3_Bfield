import numpy as np
import scipy.special as sc
import  scipy.optimize as sp 
import scipy.interpolate as interp
from scipy.interpolate import RegularGridInterpolator

def calc_params(a,zoom_zone_cu,boxsize_cu):
	'''calculates the rms velocity,density and forcing scale within cropped cube from AREPO scruct 'a' '''
	#crop cube of chosen size
	mid=boxsize_cu/2
	xmask=np.where(abs(mid-a.x)<zoom_zone_cu)	
	ymask=np.where(abs(mid-a.y)<zoom_zone_cu)
	zmask=np.where(abs(mid-a.z)<zoom_zone_cu)
	#combine masks to crop 3D cube 	
	mask=np.intersect1d(maskx,masky)
	mask=np.intersect1d(mask,maskz)
	#calculate rms velocity in cropped cube 
	v=np.sqrt(a.vx[mask]**2+a.vy[mask]**2+a.vz[mask]**2)
	v_rms=np.mean(abs(v)) * v_cu #cgs
	#mean density 
	av_rho=np.mean(a.rho) * rho_cu #cgs
	#calculate Kmin corresponding in the jeans length in cube 
	sound_speed=np.sqrt(ap.k_B.cgs.value * a.temp /(ap.m_p.cgs.value))
	jeans_length=np.sqrt(1/(ap.G.cgs.value*a.rho*rho_cu)) * sound_speed
	forcing_scale=np.mean(jeans_length) #cgs
	return v_rms,av_rho,forcing_scale

def uniform_from_dynamo(turb_type,forcing_scale,rho,v_turb):
	'''estimates the saturated integrated field strength field strength:
	-used for generating the uniform field for comparison with spectrum fields
	-refer to Schober2015 eq39
	   set turb_type as either 'burgers' or 'kolmogorov'.
	   set forcing length to boxsize or average Jeans length (cgs)
	'''

	kL=2*np.pi/forcing_scale 

	if turb_type=='burgers':
		theta=1/2
		Rm_crit=2718
		kstar=588*kL

	if turb_type=='kolmogorov':
		theta=1/3
		Rm_crit=107
		kstar=101*kL	

	gamma = (304*theta + 163) / 60
	m=3/4*gamma * (gamma*Rm_crit)**(-2*theta/(theta+1)) * (rho*v_turb**2)
	b=np.sqrt(m)
	return b


def spectrum(turb_type,forcing_scale,boxsize,N,rho,v_turb):
	'''turbulent dynamo saturated B field energy spectrum from Schober(2015).
	Sets k (2pi/lambda) limits as lambda=forcing_scale to lambda=boxsize/(N-1), where N is the 
	number of pixels you want in the real space image.
	K space will have N/2 k modes.
	Choose either 'burgers' or 'kolmogorov' for turbulence type.
	Easy option: set forcing scale as the boxsize (cgs) or use calc_params to find average Jeans length (cgs).'''
	
	kL=2*np.pi/forcing_scale

	kmin=kL
	kmax=2*np.pi / (boxsize/(N/2))
	
	if turb_type=='burgers':
		theta=1/2
		Rm_crit=2718
		kstar=588*kL
	if turb_type=='kolmogorov':
		theta=1/3
		Rm_crit=107
		kstar=101*kL
	
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


def mirror_spectrum(k,M,N):
	'''the input for inverse fourier transform requires negative k modes which come after the positive ones, 
	without repeating k=0 or k= -N/2, such that k= 0 to N/2 becomes k= 0 to N/2 straight to -(N/2 -1) to -1. 
	The spectrum M for the negative modes must mirror the postives modes
	The k array isn't inputted into the ifft function so can remain linear from 0 to (N-1) for plotting simplicity '''
	M_function=interp.interp1d(k,M,bounds_error=False,fill_value=0)
	n=int(N/2) #Niquist frequency
	k=np.linspace(0,n,n+1) #0th k mode and modes 1 to Niquist 
	k_negative=np.linspace(n+1,n+1 + n-2, n-1) #the -2 is because we dont need an extra 0th or N/2 mode 
	k_new=np.append(k,k_negative)

	M=M_function(k)
	M_reverse=M[::-1][1:-1] #reversed, but don't want 0th or N/2th mode
	M_mirror=np.append(M,M_reverse)
	return M_mirror,k_new

def modes(k,M,N):
	'''generate 3D wavevectors, aplitudes and phases, removes divergent modes
	-based on Lomax(2015)
	k is number of cycles per boxsize
	N is number of pixels in the real image, giving N/2 k modes'''
	
	#first fit a spline to power spectrum
	M_function=interp.interp1d(k,M,bounds_error=False,fill_value=0)
	
	#3D wavevectors have a 3D amplitude for each k=(kx,ky,kz)	
	#space for B field amplitude in k space
	n=int(N/2) #k modes up to Niquist frequency
	Bx=np.zeros((N,N,N),dtype=np.complex128)
	By=np.zeros((N,N,N),dtype=np.complex128)
	Bz=np.zeros((N,N,N),dtype=np.complex128)

        #values for navigating k space
	kz=np.linspace(0,N-1,N)
	ky,kx=np.meshgrid(kz,kz) #working in [x,y] format 

	print('creating k space')
	#for each coordinate in 3D k space, assign a 3D amplitude and phase
	#split k space sube into xy sheet stack
	for z in range(N): 
	
		if z==int(n):
			print('50%')
		
		#figure out where on the power spectrum we are 
		K=np.sqrt( kx**2 + ky**2 + (kz[z])**2)
		if K>0:
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
			#first normalise wavevectors
			mag=np.sqrt(kx**2+ky**2+kz[z]**2) 
			dotx=(Ax*(kx/mag))
			doty=(Ay*(ky/mag))
			dotz=(Az*(kz[z]/mag))
			dot_product=dotx+doty+dotz
			Ax=Ax-(kx/mag)*dot_product
			Ay=Ay-(ky/mag)*dot_product
			Az=Az-(kz[z]/mag)*dot_product
			mask=np.where(mag==0) #preventing nans 
			Ax[mask]=0
			Ay[mask]=0
			Az[mask]=0
			mask2=np.where(K>N-1) #K cube goes beyond k spectrum 
			Ax[mask2]=0
			Ay[mask2]=0
			Az[mask2]=0		
			#split the signal into real/imaginary parts based on phase 
			Bx[:,:,z]=Ax*(np.cos(phix)+np.sin(phix)*1j)
			By[:,:,z]=Ay*(np.cos(phiy)+np.sin(phiy)*1j)
			Bz[:,:,z]=Az*(np.cos(phiz)+np.sin(phiz)*1j)
		
	print('performing transform')	
	#reverse FFT into real space
	
	bx=np.fft.ifftn(Bx,np.array([N,N,N]))
	by=np.fft.ifftn(By,np.array([N,N,N]))
	bz=np.fft.ifftn(Bz,np.array([N,N,N]))


	return Bx,By,Bz,bx.real,by.real,bz.real



def create_field(a,turb_type,zoom_zone_cu,boxsize_cu,N):
	'''calculates turbulent v and cloud rho within cropped zoomzone, generates 
	magnetic energy density spectrum and performs 3D inverse FT to create field
	with dimentions NxNxN containing (N/2)**3 3D k modes'''
	v,rho,L=calc_params(a,zoom_zone_cu,boxsize_cu)
	boxsize=boxsize_cu*d_cu
	M,k=spectrum(turb_type,L,boxsize,N,rho,v)
	k=discrete_modes(k,M,boxsize)
	Bx,By,B,bx,by,bz=modes(k,M,N)
	return M,k,bx,by,bz



def div_check(bx,by,bz):
	'''calculate the divergence of the box using Gauss theorem
	integral (divB dV) = integral (B.n dS)'''

	L_pix=box/N
	S_pix=L_pix**2

	F_dot_S_x0 = sum(sum(bx[0,:,:]  * -S_pix))
	F_dot_S_x1 = sum(sum(bx[-1,:,:] *  S_pix))

	F_dot_S_y0 = sum(sum(by[0,:,:]  * -S_pix))
	F_dot_S_y1 = sum(sum(by[-1,:,:] *  S_pix))

	F_dot_S_z0 = sum(sum(bz[0,:,:]  * -S_pix))
	F_dot_S_z1 = sum(sum(bz[-1,:,:] * -S_pix))
	
	integral = (F_dot_S_x0+F_dot_S_x1+F_dot_S_y0+F_dot_S_y1+F_dot_S_z0+F_dot_S_z1)
	
	divB=integral / box**3 #dividing by the volume of the box

	


def rescale(bx,by,bz,theta,kstar,Rm_crit,boxsize,N,rho,v_turb):
	b=uniform_from_dynamo(theta,kstar,Rm_crit,boxsize,N,rho,v_turb) 
	energy_density=b**2 #cgs 
	E_tot_goal=energy_density*boxsize**3

	b_mag=np.sqrt(bx**2+by**2+bz**2)
	E_density_mag=b_mag**2

	pix_vol=boxsize**3/N**3
	E_per_pix=E_density_mag*pix_vol
	E_tot=sum(sum(sum(E_per_pix)))	

	E_boost_factor=E_tot_goal/E_tot
	B_boost_factor=np.sqrt(E_boost_factor)
	bx,by,bz=bx*B_boost_factor,by*B_boost_factor,bz*B_boost_factor
	return bx,by,bz




def interpolate_field(bx,by,bz,x,y,z,box,N):
	'''interpolates bx,by,bz values fromN**3 field box onto Arepo coordinates within boxsize[cgs]'''
	side=np.linspace(0,box,N)
	fx=RegularGridInterpolator((side,side,side),bx)
	fy=RegularGridInterpolator((side,side,side),by)
	fz=RegularGridInterpolator((side,side,side),bz)
	bx=fx((x,y,z))
	by=fy((x,y,z))
	bz=fz((x,y,z))
	return bx,by,bz
	
