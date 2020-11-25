import numpy as np
import matplotlib.pyplot as plt 
import code_units
from scipy.stats import binned_statistic
import arepo_utils
import astropy.constants as ap 



def read_cube(name):
        '''returns AREPO binary cube data, coordinates arranged as data[y,x,z]'''
        A=np.fromfile(name,dtype=np.float32)
        a=np.fromfile(name,dtype=np.int32)
        n=a[0]
        N=n*n
        A=A[3:]
        cube=np.zeros((n,n,n))
        for X in range(n):
                plane=A[X*N:X*N+N].reshape(n,n)
                cube[X,:,:]=plane

        return cube


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

'''||||||||| POWER SPECTRUM ||||||||||'''

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


'''||||||||||| IMAGES ||||||||||'''

def prep_image_line(dirname, file_names):
        rho={}
        for i in range(len(file_names)):
            print('reading '+dirname+file_names[i])
            rho_=read_cube(dirname+file_names[i])
            rho[i]=rho_
        return rho

def grid_plot(dirnames,file_names,xlabels,ylabels):
        '''grid of images of dimensions [len(dirnames), len(filenames)]'''

        fig,axs=plt.subplots(int(len(dirnames)),int(len(file_names)))
        for i in range(len(dirnames)):
           rho=prep_image_line(dirnames[i],file_names)
           axs[i,0].set_ylabel(ylabels[i],fontsize=15)
           for j in range(len(file_names)):
               axs[i,j].axis("off")
               axs[i,j].imshow(np.log10( np.sum(rho[j],2) / len(rho[j][:,0,0]) * code_units.rho_cu) ) 

               axs[i,j].tick_params(axis="x", labelsize=15)
               axs[i,j].tick_params(axis="y", labelsize=15)
               if i==0:
                   axs[0,j].set_title(xlabels[j],fontsize=15)
        plt.subplots_adjust(wspace=0, hspace=0)
        return fig,axs
                   

def IMF(dirs,snap,no_bins,name):
        cs='b','g','r','cyan'
        for i in range(len(dirs)):
            a=arepo_utils.aread(dirs[i]+snap)
            if i==0:
                max_sinkmass=a.sinkmass.max()
            else: 
               if a.sinkmass.max()>max_sinkmass:
                   max_sinkmass=a.sinkmass.max()
        for i in range(len(dirs)):
            #f = open('IMF'+str(8+i)+'.txt', "x")
            a=arepo_utils.aread(dirs[i]+snap)
            N,M=np.histogram(a.sinkmass,bins=np.linspace(0,max_sinkmass,no_bins))
            N_=np.array([])
            M_=np.array([])
            for j in range(len(N)):
               N_=np.append(N_,N[j])
               N_=np.append(N_,N[j])
               M_=np.append(M_,M[j])
               M_=np.append(M_,M[j+1])
            plt.plot(M_,N_,cs[i])
         #   for j in range(len(M)):
         #       f.write(str(N[j]),str(M[j]) + '\n')
         #   f.close()




def IMF_col(dirs,snap,no_bins,name):
        fig,axs=plt.subplots(4,sharex=True)
        plt.subplots_adjust(wspace=0, hspace=0)
        axs[-1].set_xlabel(r'M [M$_{\odot}$]',fontsize=20)
        #plt.ylabel( r'N$_{sink}$',fontsize=20)
        cs='b','g','r','cyan'
        rhos=r'$\rho_{sink}$=10$^{-10}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-9}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-8}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-7}$gcm$^{-3}$'
        for i in range(len(dirs)):
            a=arepo_utils.aread(dirs[i]+snap)
            if i==0:
                max_sinkmass=a.sinkmass.max()*code_units.M_cu/ap.M_sun.cgs.value
            else:
               if a.sinkmass.max()>max_sinkmass:
                   max_sinkmass=a.sinkmass.max()*code_units.M_cu/ap.M_sun.cgs.value
        for i in range(len(dirs)):
            #f = open('IMF'+str(8+i)+'.txt', "x")
            a=arepo_utils.aread(dirs[i]+snap)
            N,M=np.histogram(a.sinkmass,bins=np.linspace(0,max_sinkmass,no_bins))
            #axs[i].bar(M,N,label=rhos[i])
            axs[i].hist(a.sinkmass*code_units.M_cu/ap.M_sun.cgs.value,bins=np.linspace(0,max_sinkmass,no_bins),color=cs[i],label=rhos[i])
            #axs[i].legend(fontsize=12,frameon=False,loc='upper right')
            axs[i].set_ylim(0,N.max()+1)
            axs[i].text(0.85,0.85,rhos[i],ha='center', va='center', transform=axs[i].transAxes,fontsize=12)
            axs[i].get_yticklabels()[0].set_visible(False)
            axs[i].set_yticks([])
        return fig,axs


