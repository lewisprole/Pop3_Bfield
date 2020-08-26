import numpy as np
import matplotlib.pyplot as plt 
plt.ion()




def spectrum(N,n):
	'''Energy spectrum exponent n (not power spectrum)'''
	k=np.linspace(0,N-1,N)
	P=k**n
	return k,P

def amplitudes(k,P):
	dk=k[1]-k[0] #should be 1 
	shell_surface=4*np.pi * k**2 
	As_squared=P/shell_surface #Pdk = A_av^2 4pi k^2 dk
	return As_squared 

def kspace(k,As_squared,remove_compressive):
	Ncube=2*len(k)
	Ax_space=np.zeros((Ncube,Ncube,Ncube))
	Ay_space=np.zeros((Ncube,Ncube,Ncube))
	Az_space=np.zeros((Ncube,Ncube,Ncube))
	ks=np.fft.fftfreq(Ncube)*Ncube
	kx,ky,kz=np.meshgrid(ks,ks,ks)
	kmag=np.sqrt(kx**2+ky**2+kz**2)


	#sneaky move to avoid massive loop - flatten arrays and use kmag as the arg for A coefficients 
	Kflat=kmag.flatten().astype(int) #'int' will make all k magnitudes round down to the lower bin edge = the correct arg (dk=1)
	Axflat=Ax_space.flatten()
	Ayflat=Ay_space.flatten()
	Azflat=Az_space.flatten()
	in_range=np.where(Kflat<As_squared.shape) #the kmags go beyond the original spectrum, only use the modes up to kmag=k.max
	Ks_in_range=Kflat[in_range]
	Axflat[in_range]=np.sqrt(1/3 * As_squared[Ks_in_range])
	Ayflat[in_range]=np.sqrt(1/3 * As_squared[Ks_in_range])
	Azflat[in_range]=np.sqrt(1/3 * As_squared[Ks_in_range])
	Ax_space=Axflat.reshape(Ncube,Ncube,Ncube) #turn back into 3D array
	Ay_space=Ayflat.reshape(Ncube,Ncube,Ncube)
	Az_space=Azflat.reshape(Ncube,Ncube,Ncube)
	#^ has been tested and works^
	

	xphi=np.random.uniform(0,1,(Ncube,Ncube,Ncube))*2*np.pi
	yphi=np.random.uniform(0,1,(Ncube,Ncube,Ncube))*2*np.pi
	zphi=np.random.uniform(0,1,(Ncube,Ncube,Ncube))*2*np.pi
	
	
	if remove_compressive==True: #see Lomax 2015 section 2.1.3 
		dotx=kx*Ax_space / kmag
		doty=ky*Ay_space / kmag
		dotz=kz*Az_space / kmag
		dot_product=(dotx+doty+dotz)
		#dot_product[np.where(kmag==0)]=0
		mask=np.where(kmag>0)
		Ax_space[mask]=Ax_space[mask]- kx[mask]*dot_product[mask] /kmag[mask]
		Ax_space[mask]=Ay_space[mask]- ky[mask]*dot_product[mask] /kmag[mask]
		Ax_space[mask]=Az_space[mask]- kz[mask]*dot_product[mask] /kmag[mask]

	vxk=Ax_space*(np.cos(xphi)+1j*np.sin(xphi))
	vyk=Ay_space*(np.cos(yphi)+1j*np.sin(yphi))
	vzk=Az_space*(np.cos(zphi)+1j*np.sin(zphi))
	
	return vxk,vyk,vzk

def real_space(N,n):
	k,P=spectrum(N,n)
	As_squared=amplitudes(k,P)
	Ax,Ay,Az=kspace(k,As_squared,True)
	bx=np.fft.ifftn(Ax)
	by=np.fft.ifftn(Ay)
	bz=np.fft.ifftn(Az)
	return bx,by,bz





