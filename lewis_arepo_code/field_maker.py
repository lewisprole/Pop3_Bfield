import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator 
import scipy.special as sc
from scipy.interpolate import interp1d
plt.ion()




def spectrum(N,n):
	'''Energy spectrum exponent n (not power spectrum)'''
	k=np.linspace(0,N-1,N)
	P=k**n
	return k,P

def magnetic_spectrum(N,forcing_scale,boxsize,turb_type):
	'''magnetic energy spectrum from Schober 2015'''
	kL=2*np.pi/forcing_scale
	kmin=kL
	kmax=2*np.pi/1 #arbitrary "small" length scale=1 to get a big k range for interpolate
	k=10**np.linspace(np.log10(kL),np.log10(kmax),1000)
	if turb_type=='burgers':
		theta=1/2
		Rm_crit=2718
		kstar=588*kL
	if turb_type=='kolmogorov':
		theta=1/3
		Rm_crit=107
		kstar=101*kL
	P=k**(3/2) * sc.kn(1,k/kstar)
	f=interp1d(k,P,fill_value='extrapolate')

	#interpolate for the chosen boxsize and N 
	k=np.linspace(1,N-1,N-1) *2*np.pi /boxsize #radians per meter (just for the calculation)
	P=f(k)
	k=np.linspace(1,N-1,N-1) #cycles per boxlength (used for the rest of the field generation)

	return k,P


def amplitudes(k,P):
	dk=k[1]-k[0] #should be 1 
	shell_volume=4*np.pi * k**2 *dk
	As_squared=P/shell_volume #Pdk = A_av^2 4pi k^2 dk
	return As_squared 

def kspace(k,As_squared,remove_compressive):
	print('beginning creating phase space')
	Ncube=2*len(k) #multiply by 2 because the fft.fftfreq function sets second half of array as a mirror
	Ax_space=np.zeros((Ncube,Ncube,Ncube))
	Ay_space=np.zeros((Ncube,Ncube,Ncube))
	Az_space=np.zeros((Ncube,Ncube,Ncube))
	ks=np.fft.fftfreq(Ncube)*Ncube #puts the k array into the weird format that fft wants
	ky,kx,kz=np.meshgrid(ks,ks,ks)
	kmag=np.sqrt(kx**2+ky**2+kz**2)

	print('divide energy between dimensions')
	#sneaky move to avoid massive loop - flatten arrays and use kmag as the arg for A coefficients (since dk=1) 
	Kflat=kmag.flatten().astype(int) #'int' will make all k magnitudes round down to the lower bin edge = the correct arg (dk=1)
	Axflat=Ax_space.flatten()
	Ayflat=Ay_space.flatten()
	Azflat=Az_space.flatten()
	in_range=np.where(Kflat<As_squared.shape) #the kmags go beyond the original spectrum, only use the modes up to kmag=k.max
	Ks_in_range=Kflat[in_range]
	Axflat[in_range]=np.sqrt(1/3 * As_squared[Kflat[in_range]])
	Ayflat[in_range]=np.sqrt(1/3 * As_squared[Kflat[in_range]])
	Azflat[in_range]=np.sqrt(1/3 * As_squared[Kflat[in_range]])
	Ax_space=Axflat.reshape(Ncube,Ncube,Ncube) #turn back into 3D array
	Ay_space=Ayflat.reshape(Ncube,Ncube,Ncube)
	Az_space=Azflat.reshape(Ncube,Ncube,Ncube)
	#^ has been tested and works^

	Ax_space[0,0,0]=0 #remove NaN from 000 coordinate - Also means 0 mean field 
	Ay_space[0,0,0]=0
	Az_space[0,0,0]=0

	
	print('create phase offsets')
	xphi=np.random.uniform(0,1,(Ncube,Ncube,Ncube))*2*np.pi
	yphi=np.random.uniform(0,1,(Ncube,Ncube,Ncube))*2*np.pi
	zphi=np.random.uniform(0,1,(Ncube,Ncube,Ncube))*2*np.pi
	
	print('remove compressive modes')
	if remove_compressive=='dot': #see Lomax 2015 section 2.1.3 
		mask=np.where(kmag>0)

		dotx=kx[mask]*Ax_space[mask] / kmag[mask]
		doty=ky[mask]*Ay_space[mask] / kmag[mask]
		dotz=kz[mask]*Az_space[mask] / kmag[mask]
		dot_product=(dotx+doty+dotz)
		
		mask=np.where(kmag>0)
		projx=dot_product * kx[mask]/kmag[mask]
		projy=dot_product * ky[mask]/kmag[mask]
		projz=dot_product * kz[mask]/kmag[mask]

		Ax_space[mask]=Ax_space[mask]- projx
		Ay_space[mask]=Ay_space[mask]- projy
		Az_space[mask]=Az_space[mask]- projz

	if remove_compressive=='cross': #wrong - don't try 
		mask=np.where(kmag>0)
		Ax_space[mask]=(ky*Az_space-kz*Ay_space)[mask] /kmag[mask]
		Ay_space[mask]=-(kx*Az_space-kz*Ax_space)[mask] /kmag[mask]
		Az_space[mask]=(kx*Ay_space-ky*Ax_space)[mask] /kmag[mask]

	if remove_compressive=='no':
		print('skipping removal')



	print('apply phase offsets')
	vxk=(Ax_space*(np.cos(xphi)+1j*np.sin(xphi)))
	vyk=(Ay_space*(np.cos(yphi)+1j*np.sin(yphi)))
	vzk=(Az_space*(np.cos(zphi)+1j*np.sin(zphi)))
	
	return vxk*Ncube**3,vyk*Ncube**3,vzk*Ncube**3 #need to *N^3 because the fft output is *N^3 for some reason  



def interpolate(bx,by,bz,x,y,z,boxsize):
	'''take grid point velocities and interpolate them so that velocities can be
	given anywhere in the cube, apply to the positions of the AREPO cells'''
	print('beginning interpolation')
	size=len(bx[:,0,0])
	side=np.linspace(0,boxsize,size)

	fx=RegularGridInterpolator((side,side,side),bx)
	Bx=fx((x,y,z))

	fy=RegularGridInterpolator((side,side,side),by)
	By=fy((x,y,z))

	fz=RegularGridInterpolator((side,side,side),bz)
	Bz=fz((x,y,z))

	return Bx,By,Bz


def real_space(N,n,boxsize):
	k,P=spectrum(N,n)
	As_squared=amplitudes(k,P)
	Ax,Ay,Az=kspace(k,As_squared,'cross')
	bx=np.fft.ifftn(Ax)
	by=np.fft.ifftn(Ay)
	bz=np.fft.ifftn(Az)
	bx,by,bz=bx.real,by.real,bz.real
	x=np.linspace(0,boxsize,2*N)
	y,x,z=np.meshgrid(x,x,x)
	return bx,by,bz,x,y,z

def rescale(bx,by,bz,B_target):
	print('beginning rescale')
	bx=bx-np.mean(bx)
	by=by-np.mean(by)
	bz=bz-np.mean(bz)
	B=np.sqrt(bx**2+by**2+bz**2)
	factor=B_target/np.mean(B)
	bx=factor*bx
	by=factor*by
	bz=factor*bz
	return bx,by,bz 
	
def prepare_ICs_Bfield(N,forcing_scale,boxsize,turb_type):
	k,P=magnetic_spectrum(N,forcing_scale,boxsize,turb_type)
	plt.figure()
	plt.loglog(k,P,color='b')
	plt.loglog(k,P,'o',color='b')
	As_squared=amplitudes(k,P)
	Ax,Ay,Az=kspace(k,As_squared,'dot')
	bx=np.fft.ifftn(Ax)
	by=np.fft.ifftn(Ay)
	bz=np.fft.ifftn(Az)
	bx,by,bz=bx.real,by.real,bz.real
	return bx,by,bz

def prepare_ICs(N,n,boxsize,xarepo,yarepo,zarepo,strength):
	k,P=spectrum(N,n)
	As_squared=amplitudes(k,P)
	Ax,Ay,Az=kspace(k,As_squared,'dot')
	print('reverse transform x')
	bx=np.fft.ifftn(Ax)
	print('y')
	by=np.fft.ifftn(Ay)
	print('z')
	bz=np.fft.ifftn(Az)
	print('keeping real')
	bx,by,bz=bx.real,by.real,bz.real
	bx,by,bz=interpolate(bx,by,bz,xarepo,yarepo,zarepo,boxsize)
	bx,by,bz=rescale(bx,by,bz,strength)
	return bx,by,bz



def create_nonscaled_field(N,n):
	k,P=spectrum(N,n)
	As_squared=amplitudes(k,P)
	Ax,Ay,Az=kspace(k,As_squared,'no')
	print('reverse transform x')
	bx=np.fft.ifftn(Ax)
	print('y')
	by=np.fft.ifftn(Ay)
	print('z')
	bz=np.fft.ifftn(Az)
	print('keeping real')
	bx,by,bz=bx.real,by.real,bz.real
	return bx,by,bz
