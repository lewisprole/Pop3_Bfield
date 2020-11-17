import numpy as np
import matplotlib.pyplot as plt 
import code_units
from scipy.stats import binned_statistic

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


def subtract_radial(vx,vy,vz,x,y,z,boxsize):
        Ncube=len(vx[:,:,0])
        Lcell=boxsize/Ncube
        c=int(Ncube/2)
        rx=x[c,c,c]-x
        ry=y[c,c,c]-y
        rz=z[c,c,c]-z
        rmag=np.sqrt(rx**2+ry**2+rz**2)
        #vmag=np.sqrt(vx**2+vy**2+vz**2)
        crossx,crossy,crossz=np.zeros_like(rx),np.zeros_like(rx),np.zeros_like(rx)
        mask=np.where(rmag>0)
        crossx[mask]=(vy*rz-vz*ry)[mask]/rmag[mask]
        crossy[mask]=-(vx*rz-vz*rx)[mask]/rmag[mask]
        crossz[mask]=(vx*ry-vy*rx)[mask]/rmag[mask]
        return crossx,crossy,crossz


def create_spectrum(velfile,boxsize,subtract):
        '''reads cube, fft, creates power spectrum, please give boxsize in cm
        subtract='yes' for radial profile subtraction '''
        print('reading')
        vx,vy,vz,x,y,z=read_3cube(velfile)
        vx,vy,vz=vx*code_units.v_cu, vy*code_units.v_cu, vz*code_units.v_cu
        x,y,z=x/x.max() * boxsize, y/y.max() * boxsize, z/z.max() * boxsize

        if subtract=='yes':
           print('subtracting')
           vx,vy,vz=subtract_radial(vx,vy,vz,x,y,z,boxsize) 

        print('fft')
        Ax=np.fft.fftn(vx)
        Ay=np.fft.fftn(vy)
        Az=np.fft.fftn(vz)
        A=np.sqrt(abs(Ax)**2+abs(Ay)**2+abs(Az)**2)

        print('exploring k space')
        Ncube=int(len(A[0,0,:]))
        k=np.fft.fftfreq(Ncube)*Ncube
        kx,ky,kz=np.meshgrid(k,k,k)
        K=np.sqrt(kx**2+ky**2+kz**2)
        print('spectra')
        print('summing energies')
        bins=np.linspace(1,int(K.max()),int(K.max()))
        av1,ks1,args=binned_statistic(K.flatten(),abs(A/Ncube**3).flatten()**2,bins=bins)
        dk=ks1[1]-ks1[0]
        energy1=av1*4*np.pi*ks1[:-1]**2 *dk

        return ks1[1:],energy1



def write_txt(velfile,boxsize,subtract,name):
       ks,energy=create_spectrum(velfile,boxsize,subtract)
       f = open(name, "x")
       for i in range(len(ks)):
           f.write(str(ks[i]) + ' ' + str(energy[i]) + ' \n')
       f.close()




def txtread(txtfile):
        
        e=[]
        k=[]
        with open(txtfile) as f:
                for line in f.readlines():
                        k.append(line.split()[0])
                        e.append(line.split()[1])
        return np.asarray(k).astype(float),np.asarray(e).astype(float) 
