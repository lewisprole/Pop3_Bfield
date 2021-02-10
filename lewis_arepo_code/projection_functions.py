import numpy as np
import matplotlib.pyplot as plt 
import code_units
from scipy.stats import binned_statistic
import arepo_utils
import astropy.constants as ap
from mpl_toolkits.axes_grid1 import ImageGrid 



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




def plot_spectrum(files,subtracted_too,labels):
        if subtracted_too=='no':
           fig,axs=plt.subplots(1)
        else:
           fig,axs=plt.subplots(2,sharex=True)
           plt.subplots_adjust(hspace=0)
        for i in range(len(files)):
            k,e=txtread(files[i])
            axs[0].loglog(k,e,label=labels[i])
            axs[0].set_xlim(2,k.max()/2)
            axs[0].tick_params(axis="y", labelsize=15,direction="in")
            if subtracted_too=='yes':
                k,e=txtread(files[i]+'_radial')
                axs[1].loglog(k,e)
                axs[1].set_xlim(2,k.max()/2)
                axs[1].set_ylabel(r'$P_{v_\theta}$',fontsize=20)
                axs[1].set_xlabel('Cycles per box length',fontsize=20)
                axs[1].text(0.2,0.1,'radial profile subtracted',ha='center', va='center', transform=axs[1].transAxes,fontsize=10)
                axs[1].tick_params(axis="y", labelsize=15,direction="in")
                axs[1].tick_params(axis="x", labelsize=15,direction="in")
            else:
                axs[0].set_xlabel('Cycles per box length',fontsize=20)
                axs[0].tick_params(axis="x", labelsize=15,direction="in")
        axs[0].set_ylabel(r'$P_v$',fontsize=20)
        axs[0].legend(loc='upper right',fontsize=10,frameon=False)
        return fig,axs
            
        
        

'''||||||||||| IMAGES ||||||||||'''

def prep_image_line(dirname, file_numbers):
        rho={}
        for i in range(len(file_numbers)):
            print('reading '+dirname+'density_grid_'+file_numbers[i])
            rho_=read_cube(dirname+'density_grid_'+file_numbers[i])
            rho[i]=rho_
        return rho

def convert_sinkcoord_image(dirname,snap_no,imsize,imx1,imx2,imy1,imz1):
        '''function to convert the AREPO coordinates of sinks into the grid projection coordinates
        imsize is number of pixels per grid dimension
        imx1 and imx2 are the AREPO coords corresponding to the edges of the grids x dimension
        assumes the grid is a N**3 cube'''
        a=arepo_utils.aread(dirname+'snapshot_'+snap_no)
        X,Y,Z=a.sinkx-imx1,a.sinky-imy1,a.sinkz-imz1
        pixel=(imx2-imx1)/imsize
        x,y,z=(X/pixel).astype(int),(Y/pixel).astype(int),(Z/pixel).astype(int)
        f=open(dirname+'sink_coord_'+snap_no,'x')
        for i in range(len(np.asarray(x))):
             f.write(str(x[i])+' '+str(y[i])+' '+str(z[i])+' '+str(a.sinkmass[i])+'\n')
        f.close()
        
def read_sink(dirname,sink_number):
        x=[]
        y=[]
        z=[]
        M=[]
        with open(dirname+'sink_coord_'+sink_number) as f:
            for line in f.readlines():
                    x.append(float(line.split()[0]))
                    y.append(float(line.split()[1]))
                    z.append(float(line.split()[2]))
                    M.append(float(line.split()[3]))
        
        return np.asarray(x),np.asarray(y),np.asarray(z),np.asarray(M)

def CoM(x,y,z,M):
        X=sum(x*M)/sum(M)
        Y=sum(y*M)/sum(M)
        Z=sum(z*M)/sum(M)
        return X,Y,Z

def grid_plot(dirnames,file_numbers,xlabels,ylabels):
        '''grid of images of dimensions [len(dirnames), len(filenames)]'''

        fig,axs=plt.subplots(int(len(dirnames)),int(len(file_numbers)))

        for i in range(len(dirnames)):
           rho=prep_image_line(dirnames[i],file_numbers)
           for j in range(len(file_numbers)):
               axs[i,j].axis("off")
               axs[i,j].imshow(np.log10( np.sum(rho[j],2) / len(rho[j][:,0,0]) * code_units.rho_cu) )

               mid=np.where(rho[j]==rho[j].max())
               if len(mid[0])>1:
                   axs[i,j].set_ylim(mid[0][0]-200,mid[0][0]+200)
                   axs[i,j].set_xlim(mid[1][0]-200,mid[1][0]+200)
               else:
                   axs[i,j].set_ylim(mid[0]-200,mid[0]+200)
                   axs[i,j].set_xlim(mid[1]-200,mid[1]+200)


               axs[i,j].tick_params(axis="x", labelsize=15)
               axs[i,j].tick_params(axis="y", labelsize=15)
               if i==0:
                   axs[0,j].set_title(xlabels[j],fontsize=15)
               if j==0:
                   axs[i,0].set_ylabel(ylabels[i],fontsize=15)
        plt.subplots_adjust(wspace=-0.6, hspace=0)
        return fig,axs


def grid_with_sinks(dirnames,file_numbers,xlabels,ylabels):
        fig = plt.figure(figsize=(len(dirnames), len(file_numbers)))
        grid = ImageGrid(fig, 111,  # similar to subplot(111)
		nrows_ncols=(5, 4),  # creates 2x2 grid of axes
		axes_pad=0,
		cbar_mode="single",  # pad between axes in inch.
		 ) 
        for i in range(len(dirnames)):
           I=int(len(dirnames) - i -1)
           rho=prep_image_line(dirnames[I],file_numbers)
           for j in range(len(file_numbers)):
               J=int(len(file_numbers) - j -1)
               print(I*len(file_numbers)+J)
               grid[I*len(file_numbers)+J].set_yticks([])
               grid[I*len(file_numbers)+J].set_xticks([])
               if i==0:
                    if j==0:
                        vmin=(np.log10( np.sum(rho[J],2) / len(rho[J][:,0,0]) * code_units.rho_cu)).min()
                        vmax=(np.log10( np.sum(rho[J],2) / len(rho[J][:,0,0]) * code_units.rho_cu)).max()
               #import the sinks
               x,y,z,M=read_sink(dirnames[I],file_numbers[J])
               mask=np.where((np.sqrt((y-500)**2+ (x-500)**2) <400))
               cx,cy,cz=CoM(x[mask],y[mask],z[mask],M[mask])

               data=np.log10( np.sum(rho[J],2) / len(rho[J][:,0,0]) * code_units.rho_cu)
               im=grid[I*len(file_numbers)+J].imshow(data[int(cx)-200:int(cx)+200,int(cy)-200:int(cy)+200],cmap='bone',vmin=vmin,vmax=vmax)

               mask=np.where((np.sqrt((y-cy)**2+ (x-cx)**2) <200))
               grid[I*len(file_numbers)+J].scatter(y[mask]-(cy-200),x[mask]-(cx-200),s=0.5,c='magenta')
               grid[I*len(file_numbers)+J].text(0.75,0.08,r'N$_{sink}$='+str(len(x)),ha='center', va='center', transform=grid[I*len(file_numbers)+J].transAxes,fontsize=5,color='w')

               grid[I*len(file_numbers)+J].tick_params(axis="x", labelsize=5)
               grid[I*len(file_numbers)+J].tick_params(axis="y", labelsize=5)
               if I==0:
                   grid[I*len(file_numbers)+J].set_title(xlabels[J],fontsize=5)
               if J==0:
                   grid[I*len(file_numbers)+J].set_ylabel(ylabels[I],fontsize=5,rotation=75)
        cbar=grid.cbar_axes[0].colorbar(im)
        cbar.ax.tick_params(labelsize=5)
        cbar.ax.set_ylabel(r'Log$_{10}$($\rho$ [gcm$^{-3}$])', rotation=270,fontsize=5,labelpad=15)
        #plt.subplots_adjust(wspace=-0.735, hspace=0)
        return fig,axs
        



def duel_column(dirs,nos):
        '''graph for 2 rows of 3 panels, showing 3 resolution criteria performed for 2 sink creation densities'''
        fig = plt.figure(figsize=(3,2))
        grid = ImageGrid(fig, 111,  # similar to subplot(111)
                nrows_ncols=(3, 2),  # creates 2x2 grid of axes
                axes_pad=(0.2,0),
                cbar_mode="single",  
                cbar_pad=0
                 )
        xlabels='8 cells','16 cells','32 cells'
        grid[0].set_title(r'$\rho_{\rm sink}$=10$^{-10}$gcm$^{-3}$'
                             '\n'
                              r't=2500yrs',fontsize=10)
        grid[1].set_title(r'$\rho_{\rm sink}$=10$^{-7}$gcm$^{-3}$'
                             '\n'
                              r't=350yrs',fontsize=10)
        for i in range(len(nos)):
                
                grid[i*2].set_ylabel(xlabels[i],fontsize=10)

                rho=read_cube(dirs[1]+'density_grid_'+nos[i])
                data=np.log10( np.sum(rho,2) / len(rho[:,0,0]) * code_units.rho_cu)
                x,y,z,M=read_sink(dirs[1],nos[i])
                mask=np.where((np.sqrt((y-500)**2+ (x-500)**2) <400))
                cx,cy,cz=CoM(x[mask],y[mask],z[mask],M[mask])
                cx,cy,cz=int(cx),int(cy),int(cz)
                if i==0:
                       cx,cy,cz=600,600,600
                vmin,vmax=data.min(),data.max()
                im= grid[i*2+1].imshow(data[cx-400:cx+400,cy-400:cy+400],vmin=vmin,vmax=vmax,cmap='bone')
                mask=np.where((np.sqrt((y-cy)**2+ (x-cx)**2) <400))
                grid[i*2+1].scatter(y[mask]-(cy-400),x[mask]-(cx-400),s=0.5,c='magenta')
                grid[i*2+1].text(0.8,0.08,r'N$_{sink}$='+str(len(x)),ha='center', va='center', transform=grid[i*2+1].transAxes,fontsize=10,color='w')

                rho=read_cube(dirs[0]+'density_grid_'+nos[i])
                data=np.log10( np.sum(rho,2) / len(rho[:,0,0]) * code_units.rho_cu)
                x,y,z,M=read_sink(dirs[0],nos[i])
                mask=np.where((np.sqrt((y-500)**2+ (x-500)**2) <400))
                cx,cy,cz=CoM(x[mask],y[mask],z[mask],M[mask])
                cx,cy,cz=int(cx),int(cy),int(cz)
                if i==2:
                       cx,cy,cz=500,500,500
                grid[i*2].imshow(data[cx-400:cx+400,cy-400:cy+400],vmin=vmin,vmax=vmax,cmap='bone')
                mask=np.where((np.sqrt((y-cy)**2+ (x-cx)**2) <400))
                grid[i*2].scatter(y[mask]-(cy-400),x[mask]-(cx-400),s=0.5,c='magenta')
                grid[i*2].set_xlim(0,800)
                grid[i*2].set_ylim(800,0)
                grid[i*2].text(0.8,0.08,r'N$_{sink}$='+str(len(x)),ha='center', va='center', transform=grid[i*2].transAxes,fontsize=10,color='w')

                grid[i*2].tick_params(axis="x", labelsize=9)
                grid[i*2].tick_params(axis="y", labelsize=9)
                grid[i*2].set_yticks([])
                grid[i*2].set_xticks([])

                grid[i].tick_params(axis="x", labelsize=9)
                grid[i].tick_params(axis="y", labelsize=9)
                grid[i].set_yticks([])
                grid[i].set_xticks([])


                cbar=grid.cbar_axes[0].colorbar(im)
                cbar.ax.tick_params(labelsize=9)
                cbar.ax.set_ylabel(r'Log$_{10}$($\rho$ [gcm$^{-3}$])', rotation=270,fontsize=10,labelpad=20)



def IMF_col(dirs,snap,no_bins):
        fig,axs=plt.subplots(5,sharex=True)
        plt.subplots_adjust(wspace=0, hspace=0)
        axs[-1].set_xlabel(r'M [M$_{\odot}$]',fontsize=11)
        axs[-1].tick_params(axis="x", labelsize=10)
        axs[2].set_ylabel(r'N$_{\rm M}    $',fontsize=11,rotation=0)
        axs[2].yaxis.set_label_coords(-0.15,0.3)
        axs[1].yaxis.set_label_coords(-0.08, 0.1)
        #plt.ylabel( r'N$_{sink}$',fontsize=20)
        cs='b','g','r','cyan','Purple'
        rhos=r'$\rho_{sink}$=10$^{-10}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-9}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-8}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-7}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-6}$gcm$^{-3}$'
        for i in range(len(dirs)):
            a=arepo_utils.aread(dirs[i]+str(snap[i]))
            if i==0:
                max_sinkmass=a.sinkmass.max()*code_units.M_cu/ap.M_sun.cgs.value
                #min_sinkmass=a.sinkmass.min()*code_units.M_cu/ap.M_sun.cgs.value
            else:
               if a.sinkmass.max()>max_sinkmass:
                   max_sinkmass=a.sinkmass.max()*code_units.M_cu/ap.M_sun.cgs.value
               #if a.sinkmass.min()<min_sinkmass:
                   #min_sinkmass=a.sinkmass.min()*code_units.M_cu/ap.M_sun.cgs.value

        min_sinkmass=4/3*np.pi*(1.71E-05*code_units.d_cu)**3 *1e12*code_units.rho_cu /ap.M_sun.cgs.value
        bins=10**np.linspace(np.log10(min_sinkmass),np.log10(max_sinkmass),no_bins)
        for i in range(len(dirs)):
            #f = open('IMF'+str(8+i)+'.txt', "x")
            a=arepo_utils.aread(dirs[i]+str(snap[i]))
            N,M=np.histogram(a.sinkmass*code_units.M_cu/ap.M_sun.cgs.value,bins=bins)
            #axs[i].bar(M,N,label=rhos[i])
            axs[i].hist(a.sinkmass*code_units.M_cu/ap.M_sun.cgs.value,bins=bins,color=cs[i],label=rhos[i])
            #axs[i].legend(fontsize=12,frameon=False,loc='upper right')
            axs[i].set_ylim(0,N.max()+1)
            minmass=np.array([1e8,1e9,1e10,1e11,1e12])[i] *code_units.rho_cu  * 4/3*np.pi * (np.array([0.001376823,0.0004563,0.000152667,5.04801E-05,1.71E-05])[i]*code_units.d_cu)**3 /ap.M_sun.cgs.value
            axs[i].axvline(x=minmass,ymin=0,ymax=1,color='k',linestyle='--')
            axs[i].tick_params(axis="y", labelsize=10)
            axs[i].set_ylim(0,N.max()+1)
            axs[i].text(1.22,0.5,rhos[i],ha='center', va='center', transform=axs[i].transAxes,fontsize=10)
            #axs[i].get_yticklabels()[0].set_visible(False)
            axs[i].set_yticks([N.max()])
            #axs[i].set_yticks([])
            axs[i].set_xscale('log')
            plt.subplots_adjust(left = 0.15,bottom = 0.17,right=0.7)
        return fig,axs

