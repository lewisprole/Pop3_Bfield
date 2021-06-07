import numpy as np
import field_maker 
from scipy.stats import binned_statistic
import matplotlib.pyplot as plt

'''couple of functions for reading the field file and then create 1D power spectrum'''

def read_3cube(name):
        '''returns 3 cubes giving x,y and z components, individual cubes have positional
        axis cube[y,x,z]'''
        A=np.fromfile(name,dtype=np.float32)
        a=np.fromfile(name,dtype=np.int32)
        n=a[0]
        N=n*n
        A=A[3:]
        Ax=A[0::3]
        Ay=A[1::3]
        Az=A[2::3]
        cubex=np.zeros((n,n,n))
        cubey=np.zeros((n,n,n))
        cubez=np.zeros((n,n,n))
        for X in range(n):
                planex=Ax[X*N:X*N+N].reshape(n,n)
                cubex[X,:,:]=planex
                planey=Ay[X*N:X*N+N].reshape(n,n)
                cubey[X,:,:]=planey
                planez=Az[X*N:X*N+N].reshape(n,n)
                cubez[X,:,:]=planez
        x=np.linspace(0,n-1,n)
        y,x,z=np.meshgrid(x,x,x)
        return cubex,cubey,cubez,x,y,z

def create_spectrum(vx,vy,vz):
        '''reads cube, fft, creates power spectrum, please give boxsize in cm
        subtract='yes' for radial profile subtraction '''
        print('fft')
        Ax=np.fft.fftn(vx)
        Ay=np.fft.fftn(vy)
        Az=np.fft.fftn(vz)
        A=np.sqrt(abs(Ax)**2+abs(Ay)**2+abs(Az)**2)
        print('exploring k space')
        Ncube=int(len(A[0,0,:]))
        k=np.fft.fftfreq(Ncube)*(Ncube)
        kx,ky,kz=np.meshgrid(k,k,k)
        K=np.sqrt(kx**2+ky**2+kz**2)
        print('spectra')
        print('summing energies')
        bins=np.linspace(1,int(K.max()),int(K.max()))
        av1,ks1,args=binned_statistic(K.flatten(),abs(A/Ncube**3).flatten()**2,bins=bins)
        dk=ks1[1]-ks1[0]
        energy1=av1*4*np.pi*ks1[:-1]**2 *dk

        #get rid of the artifacts 
        ks1=ks1[:int(Ncube/2)+1]
        energy1=energy1[:int(Ncube/2)+1]
        return ks1,energy1

def interp_graph(files):
	labels=r'0.05N$_{\rm field}$',r'0.25N$_{\rm field}$',r'0.5N$_{\rm field}$',r'N$_{\rm field}$',r'2N$_{\rm field}$',r'10N$_{\rm field}$',r'100N$_{\rm field}$'
	#first make spectrum of original uniform field
	#bx,by,bz=field_maker.create_nonscaled_Bfield(50,3/2)
	#k,P=create_spectrum(bx,by,bz)
	#plt.loglog(k,P/k**(3/2),label=r'k$^{3/2}$',color='k',linestyle='--')
	plt.axhline(y=1.66e-6,color='k',linestyle='--')
	plt.axhline(y=1.06e-6,color='k',linestyle='--')
	for i in range(len(files)):
		bx,by,bz,x,y,z=read_3cube(files[i])
		k,P=create_spectrum(bx,by,bz)
		plt.loglog(k,P/k**(3/2),label=labels[i])
	plt.legend(fontsize=8,frameon=False)
	plt.ylabel(r'P/$k^{3/2}$',rotation=0)
	plt.xlabel('k')
	plt.ylim(5.5e-7,1.5e-5)
	plt.xlim(0.995,25.1)
	plt.tick_params(axis="x", labelsize=10,direction="in",which='both')
	plt.tick_params(axis="y", labelsize=10,direction="in",which='both')
	
