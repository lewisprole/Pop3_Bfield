import numpy as np
import matplotlib.pyplot as plt
import code_units
import sys
import io
import arepo_utils
import astropy.constants as ap
from  scipy import interpolate
import analyse
from scipy.stats import binned_statistic
plt.ion()

def Rs(a):
	mid=np.where(a.rho==a.rho.max())
	rs=np.sqrt((a.x-a.x[mid])**2+(a.y-a.y[mid])**2+(a.z-a.z[mid])**2)
	return rs


def solenoidal_velocity(a):
        '''calculates radial and non radial component of the velocity, also returns total vel magnitude'''
        #create vector to center
        mid=np.where(a.rho==a.rho.max())
        xmid=a.x[mid]
        ymid=a.y[mid]
        zmid=a.z[mid]
        r=np.array([a.x-xmid, a.y-ymid, a.z-zmid])
        #radial and non-radial velocities
        dot_product=(r[0]*a.vx + r[1]*a.vy + r[2]*a.vz)
        cross_product=(   (r[1]*a.vz - r[2]*a.vy) , - (r[0]*a.vz - r[2]*a.vx) , (r[0]*a.vy - r[1]*a.vx)   )

        #normalise
        R=np.sqrt(r[0]**2+r[1]**2+r[2]**2)
        radial_v=abs(dot_product/R)
        nonrad_v=abs(cross_product/R)
        #remove NaNs
        radial_v[np.where(R==0)]=0
        nonrad_v[:,np.where(R==0)]=np.zeros((3,1,1))
        #reduce cross product to 1D
        nonrad_v=np.sqrt(nonrad_v[0]**2 + nonrad_v[1]**2 + nonrad_v[2]**2)

        v_mag=np.sqrt(a.vx**2+a.vy**2+a.vz**2)
        return radial_v,nonrad_v,v_mag,R

def abundances(a):
	H2=a.chem[:,0]
	e=a.chem[:,1]+a.chem[:,2]+a.chem[:,4]+2*a.chem[:,5]
	return H2,e
	

def weighter(a,type):

	if type=='volume':
		normaliser=a.mass/a.rho
	if type=='mass':
		normaliser=a.mass
	return normaliser


def radial_average(variable,a,weight_type,bins,cumulative):
	rs=Rs(a)
	x=np.linspace(np.log10(np.sort(rs)[1]),np.log10(rs.max()),bins)
	x=10**x
	weight=weighter(a,weight_type)
	avs=np.array([])
	for i in range(bins-1):
		if cumulative==True:
			mask=np.where((rs>x[0]) & (rs<x[i+1]))
		else:
			mask=np.where((rs>x[i]) & (rs<x[i+1]))
		RS=rs[mask]
		if sum(weight[mask])>0:
			av=sum(variable[mask]*weight[mask])/sum(weight[mask])
		else:
			av=0
		avs=np.append(avs,av)
	return x,avs
	
def cycle_plot(snaps,weight_type,bins,labels,title):
	
	fig,axs=plt.subplots(6,sharex=True)
	plt.subplots_adjust(wspace=0, hspace=0,top=0.95,bottom=0.08,right=0.97,left=0.2)
	for i in range(len(snaps)):
		a=arepo_utils.aread(snaps[i])
		rv,nv,v,rs=solenoidal_velocity(a)
		H2,e=abundances(a)

		print('averaging')
		x,v=radial_average(v,a,weight_type,bins,False)
		print('v')
		x,rho=radial_average(a.rho,a,weight_type,bins,False)
		print('rho')
		x,rv=radial_average(rv,a,weight_type,bins,False)
		print('radial v')
		x,nv=radial_average(nv,a,weight_type,bins,False)
		print('rotational v')
		x,T=radial_average(a.temp,a,weight_type,bins,False)
		print('T')
		x,H2=radial_average(H2,a,weight_type,bins,False)
		print('H2')
		x,e=radial_average(e,a,weight_type,bins,False)
		print('e')
		

		rho=rho*code_units.rho_cu
		v=v*code_units.v_cu/1e5
		rv=rv*code_units.v_cu/1e5
		nv-nv*code_units.v_cu/1e5
		x=x*code_units.d_cu/ap.pc.cgs.value
		
		x=x[:-1]
		axs[0].loglog(x,rho,label=labels[i])
		axs[1].semilogx(x,rv)
		axs[2].semilogx(x,nv)
		axs[3].loglog(x,T)
		axs[4].loglog(x,H2)
		axs[5].loglog(x,e)
	
	axs[5].set_xlabel('R [pc]',fontsize=11)
	axs[5].tick_params(axis="x", labelsize=11,direction="in")

	axs[0].legend(fontsize=9,frameon=False,loc='lower left',bbox_to_anchor=(0.,-0.05))	
	axs[0].set_ylabel(r'$\rho$ [gcm$^{-3}$]',fontsize=11)
#	logrange=np.log10(rho.max())-np.log10(rho.min())
#	p=np.round( np.log10(rho.min()) + logrange*np.array([0.25,0.5,0.75]) ,0)
#	axs[0].set_yticks(10**(p.astype(float)))
	axs[0].set_yticks([10**-21,10**-18,10**-15,10**-12])
	axs[0].tick_params(axis="y", labelsize=11,direction="in")

	axs[1].set_ylabel(r'$v_{r}$ [kms$^{-1}$]',fontsize=11)	
#	axs[1].set_yticks((rv.min() + (rv.max()-rv.min()) * np.array([0.25,0.5,0.75])).astype(int))
	axs[1].set_yticks([1,2,3,4])
	axs[1].tick_params(axis="y", labelsize=11,direction="in")

	axs[2].set_ylabel(r'$v_{\theta}$ [kms$^{-1}$]',fontsize=11)
#	t=axs[2].get_yticks()
#	axs[2].set_yticks((nv.min() + (nv.max()-nv.min()) * np.array([0.25,0.5,0.75])).astype(int))
	axs[2].set_yticks([3,5,7,9])
	axs[2].tick_params(axis="y", labelsize=11,direction="in")


	axs[3].set_ylabel('T [K]',fontsize=11)
#	logrange=np.log10(T.max())-np.log10(T.min())
#	p=np.round( np.log10(T.min()) + logrange*np.array([0.25,0.5,0.75]) ,0)
#	axs[3].set_yticks(10**(p.astype(float)))
	axs[3].set_yticks([10**3])
	axs[3].tick_params(axis="y", labelsize=11,direction="in")
	
	axs[4].set_ylabel(r'H$_{2}$',fontsize=11)
#	logrange=np.log10(H2.max())-np.log10(H2.min())
#	p=np.round( np.log10(H2.min()) + logrange*np.array([0.25,0.5,0.75]) ,0)
#	axs[4].set_yticks(10**(p.astype(float)))
	axs[4].set_yticks([10**-2,10**-1])
	axs[4].tick_params(axis="y", labelsize=11,direction="in")

	axs[5].set_ylabel(r'$e^{-}$',fontsize=11)
#	logrange=np.log10(e.max())-np.log10(e.min())
#	p=np.round( np.log10(e.min()) + logrange*np.array([0.25,0.5,0.75]) ,0)
#	axs[5].set_yticks(10**(p.astype(float)))
	axs[5].set_yticks([10**-9,10**-8,10**-7])
	axs[5].tick_params(axis="y", labelsize=11,direction="in")
	

	for i in range(6):
		axs[i].tick_params(axis="x", labelsize=11,direction="in")
		axs[i].get_yaxis().set_label_coords(-0.15,0.5)

	plt.xlim(x.min(),np.sort(x)[-2])
	axs[0].set_title(title,fontsize=11)
	return fig,axs


def cumulative_velocity(snaps,weight_type,bins,labels):
	fig,axs=plt.subplots(3,sharex=True)
	plt.subplots_adjust(wspace=0, hspace=0,top=0.95,bottom=0.15,right=0.97,left=0.2)
	for i in range(len(snaps)):
		a=arepo_utils.aread(snaps[i])
		rv,nv,v,rs=solenoidal_velocity(a)
		print('averaging')
		x,v=radial_average(v,a,weight_type,bins,True)
		print('v')
		x,rv=radial_average(rv,a,weight_type,bins,True)
		print('rv')
		x,nv=radial_average(nv,a,weight_type,bins,True)
		print('nv')
		v,rv,nv=v*code_units.v_cu/1e5,rv*code_units.v_cu/1e5,nv*code_units.v_cu/1e5
		x=x*code_units.d_cu/ap.pc.cgs.value
		axs[0].semilogx(x[:-1],v,label=labels[i])
		axs[1].semilogx(x[:-1],rv)
		axs[2].semilogx(x[:-1],nv)
	axs[0].legend(fontsize=9,frameon=False,loc='lower center',bbox_to_anchor=(0.5,-0.05))
	axs[0].set_ylabel(r'$\bar{v}$ [kms$^{-1}$]',fontsize=11)
	axs[1].set_ylabel(r'$\bar{v_{r}}$ [kms$^{-1}$]',fontsize=11)
	axs[2].set_ylabel(r'$\bar{v_{\theta}}$ [kms$^{-1}$]',fontsize=11)
	axs[2].set_xlabel(r'R [pc]',fontsize=11)
	plt.xlim(x.min(),np.sort(x)[-2])
	for i in range(3):
		axs[i].tick_params(axis="x", labelsize=11,direction="in")
		axs[i].tick_params(axis="y", labelsize=11,direction="in")
		axs[i].get_yaxis().set_label_coords(-0.1,0.5)

	#specific to this graph
	axs[0].set_yticks([2,4,6])
	axs[1].set_yticks([1,2,3,4])
	axs[2].set_yticks([1,2,3])
	return fig,axs

def cumulative_radial(variable,a,bins):
	rs=Rs(a)
	x=np.linspace(np.log10(np.sort(rs)[1]),np.log10(rs.max()),bins)
	x=10**x
	tots=[]
	for i in range(bins-1):
		mask=np.where(rs<x[i+1])
		tots.append(sum(variable[mask]))
	return x,tots

def total_energy(snaps,weight_type,bins,labels):
	fig,axs=plt.subplots(3,sharex=True)
	plt.subplots_adjust(wspace=0, hspace=0,top=0.95,bottom=0.15,right=0.97,left=0.2)
	for i in range(len(snaps)):
		a=arepo_utils.aread(snaps[i])
		rv,nv,v,rs=solenoidal_velocity(a)
		Er,En,E=0.5*(a.mass*code_units.M_cu)*(rv*code_units.v_cu)**2, 0.5*(a.mass*code_units.M_cu)*(nv*code_units.v_cu)**2, 0.5*(a.mass*code_units.M_cu)*(v*code_units.v_cu)**2
		print('averaging')
		x,E=cumulative_radial(E,a,bins)
		print('v')
		x,Er=cumulative_radial(Er,a,bins)
		print('rv')
		x,En=cumulative_radial(En,a,bins)
		print('nv')
		
		axs[0].semilogx(x[:-1],E,label=labels[i])	
		axs[1].semilogx(x[:-1],Er)
		axs[2].semilogx(x[:-1],En)
	
	axs[0].legend(fontsize=9,frameon=False)
	axs[0].set_ylabel(r'$E_{tot}$ [ergs]',fontsize=11)
	axs[1].set_ylabel(r'$E_{r}$ [ergs]',fontsize=11)
	axs[2].set_ylabel(r'$E_{theta}$ [ergs]',fontsize=11)
	axs[2].set_xlabel(r'R [pc]',fontsize=11)
	return fig,axs

'''		|||||||||||||||||||||||||||||||| Functions for handling Fourier transforms of the velocity fields (OBSOLETE!) ||||||||||||||||||||||||||||||||||		'''

def FT_reduce1D(v_of_k,bins):
	'''reduces 3D k space into 3D spectrum'''
	x=np.linspace(0,len(v_of_k[:,0,0]),len(v_of_k[:,0,0])-1)
	y,x,z=np.meshgrid(x,x,x)
	K=np.sqrt(x**2+y**2+z**2)
	ks=np.linspace(1,x.max(),bins)
	
	vs=np.array([])
	for i in range(bins-1):
		mask=np.where((K>ks[i]) & (K<ks[i+1]))
		vs=np.append(vs,np.mean(v_of_k[mask]))
	return ks[:-1],vs

def spectrum_fit(ks,vs,K):
	'''interpolate to extract v(k) for given k'''
	f=interpolate.interp1d(ks,vs)
	return f(K)

def convert_L_to_k(scale,boxsize):
	'''convert length sclae into cycles per boxsize'''
	k= 2*np.pi/scale *boxsize/(2*np.pi) 
	return k

def jeans_length(rho,T):
	'''calculate jeans length'''
	sound_speed=np.sqrt(ap.k_B.cgs.value * T /(ap.m_p.cgs.value))
	jeans_length=np.sqrt(np.pi / ap.G.cgs.value / (rho*code_units.rho_cu)) * sound_speed /code_units.d_cu
	return jeans_length

def crop(a,croplength):
	mid=np.where(a.rho==a.rho.max())
	maskx=np.where((a.x > a.x[mid]-croplength) & (a.x > a.x[mid]+croplength))
	masky=np.where((a.y > a.y[mid]-croplength) & (a.y > a.y[mid]+croplength))
	maskz=np.where((a.z > a.z[mid]-croplength) & (a.z > a.z[mid]+croplength))
	MASK=np.intersect1d(maskx,masky)
	MASK=np.intersect1d(MASK,maskz)
	return MASK

def jeans_v(cubefile,snapshot,boxsize,bins,labels):
	fig,axs=plt.subplots(1)
	#different crop sizes
#	rs=np.array([1e-3,1e-2,1e-1]) *ap.pc.cgs.value/code_units.d_cu #code units
	for i in range(len(cubefile)):
#		rs=1e-2*ap.pc.cgs.value/code_units.d_cu
		vx,vy,vz=analyse.read_3cube(cubefile[i])
		Ncube=len(vx[:,:,0])
		a=arepo_utils.aread(snapshot[i])
		jeans=1e-4*ap.pc.cgs.value/code_units.d_cu
#	jeans=jeans_length(a.rho,a.temp)
#	JEANSV=[]
	
#		print('starting crop '+str(rs))
#		mask=crop(a,rs)
#		jeans_k=convert_L_to_k(jeans,2*rs)

#		Lcell=boxsize/Ncube		
#		Ncell=int(np.round(rs/Lcell,0))
#		Nhalf=int(np.round(Ncube/2,0))
		#amplitudes from FT of cropped v cube
		print('starting transform')
		A=np.fft.fftn(vx)#[Nhalf-Ncell:Nhalf+Ncell+1, Nhalf-Ncell:Nhalf+Ncell+1, Nhalf-Ncell:Nhalf+Ncell+1])
		normalise=int((Ncube)**3)#(2*Ncell)**2)
		print('reducing to 1D')
		ks,vs=FT_reduce1D(abs(A)/normalise,bins)
#		jeans_v=spectrum_fit(ks,vs,jeans_k)


#		axs.loglog(2*np.pi/(ks*2*np.pi/(2*Ncell*Lcell)) *code_units.d_cu/ap.pc.cgs.value ,vs*code_units.v_cu/1e5,label=labels[i])
		axs.loglog(2*np.pi/(ks*2*np.pi/(boxsize)) *code_units.d_cu/ap.pc.cgs.value ,vs*code_units.v_cu/1e5,label=labels[i])
	axs.set_xlabel(r'$\lambda$ [pc]',fontsize=11)
	axs.set_ylabel(r'$v(\lambda)$ [kms$^{-1}$]',fontsize=11)
	axs.legend(fontsize=9,frameon=False)
#		JEANSV.append(jeans_v)
	
	return fig,axs



'''             |||||||||||||||||||||||||||||||||||| Functions for handling Fourier transforms of the velocity fields NEW ||||||||||||||||||||||||||||||||||||||                  '''

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


def power_spectrum(A,boxsize,interval):
	powers=np.array([]) #space for power spectrum
	Ncube=len(A[:,:,0]) #length of side of k space cube
	x=np.linspace(0,Ncube-1,Ncube) #k modes
	k=x * 2*np.pi / boxsize #convert k modes into k=2pi/lambda
	dk=k[1]-k[0] #will be needed when integrating the spectrum
	ky,kx,kz=np.meshgrid(k,k,k) #k coordinates for navigating in k-space
	
	k_walk=k[1::interval]#np.linspace(k[1],k[-1],bins) #going to integrate outward along the 3D cube (excluding 0 - ruins log)
	
	k_mag=np.sqrt(kx**2+ky**2+kz**2) #each k mode has a radius in k space from 0
	for i in range(int(len(k_walk)-1)):
		print(str(i) +' of ' +str(int(len(k)/interval)))
		delta_k=k_walk[i+1]-k_walk[i]
		mask=np.where((k_mag>=k_walk[i]) & (k_mag<k_walk[i+1])) #shell of thickness 
		power=sum(A[mask]**2*dk**3) /delta_k #power = energy density * volume, treat v**2 as energy density in k space
					#integrate the shell by summing each 3d k mode with its dk^3 
		if delta_k<dk:
			print('dk: '+str(dk)+', delta_k: '+str(delta_k) +': '+str(k_walk[i+1]) +'-' + str(k_walk[i]))
	#		powers=np.append(powers,np.nan)			#only integrating by dk^2 becaue want the final units to be P(v)dk=E
	#	else:
		powers=np.append(powers,power)
	return k_walk[:-1],powers

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
	
def create_spectrum(A):
	Ncube=int(len(A[0,0,:]))
	k=np.fft.fftfreq(Ncube)*Ncube
	kx,ky,kz=np.meshgrid(k,k,k)
	K=np.sqrt(kx**2+ky**2+kz**2)
	bins=np.linspace(1,int(K.max()),int(K.max())+1)
	av,ks,args=binned_statistic(K.flatten(),(abs(A)/Ncube**3).flatten()**2,bins=bins)
	dk=ks[1]-ks[0]
	energy=av*4*np.pi*ks[1:]**2 *dk
	return ks[1:],energy
		
	
def cycle_spectrum(cubefiles,boxsize,interval,labels):
	'''give boxsize in cm'''
	sizes=np.array([8,6,4,2])
	marker=np.array(['v','.','+','2'])
	fig,axs=plt.subplots(2,sharex=True)
	for i in range(len(cubefiles)):	
		print('reading')
		vx1,vy1,vz1,x,y,z=read_3cube(cubefiles[i])
		vx1,vy1,vz1=vx1*code_units.v_cu, vy1*code_units.v_cu, vz1*code_units.v_cu
		x,y,z=x/x.max() * boxsize, y/y.max() * boxsize, z/z.max() * boxsize
		for j in range(2):
			if j==1:
				print('subtracting radial')
				vx,vy,vz=subtract_radial(vx1,vy1,vz1,x,y,z,boxsize)
			else:
				vx,vy,vz=vx1,vy1,vz1
			Ncube=len(vx[:,:,0])
			print('vx transform')
			Ax=np.fft.fftn(vx)
			Ax=Ax[:int(Ncube/2),:int(Ncube/2),:int(Ncube/2)]
			print('vy transform')
			Ay=np.fft.fftn(vy)
			Ay=Ay[:int(Ncube/2),:int(Ncube/2),:int(Ncube/2)]
			print('vz transform')
			Az=np.fft.fftn(vz)
			Az=Az[:int(Ncube/2),:int(Ncube/2),:int(Ncube/2)]

			A=np.sqrt(abs(Ax)**2+abs(Ay)**2+abs(Az)**2)
			#normalise=int((Ncube)**3)
			#print('normalising')
			#A=abs(A)/normalise
			print('creating power spectrum')
			#k,P=power_spectrum(A,boxsize,interval)
			k,P=create_spectrum(A)
			if j==0:
				axs[j].loglog(k,P,marker[i],markersize=4,label=labels[i])
			else:
				axs[j].loglog(k,P,marker[i],markersize=4)
			axs[j].set_ylabel(r'$P_{v}$',fontsize=11)
			axs[j].tick_params(axis="y", labelsize=11,direction="in")
	axs[1].tick_params(axis="x", labelsize=11,direction="in")
	axs[1].set_xlabel(r'$cycles per box length$]',fontsize=11)
	#axs[1].annotate('Radial profile subtracted',(3e-19,0.5e-33))
	#axs[1].loglog(k[np.where((k>2e-17) & (k<1e-16))],k[np.where((k>2e-17)&(k<1e-16))]**(-5/3)*2e-60,linestyle='--',color='k')
	#axs[1].annotate(r'$\propto k^{-5/3}$',(4e-17,1e-32))
	#axs[1].loglog(k[np.where((k>2e-18) & (k<2e-17))],k[np.where((k>2e-18) & (k<2e-17))]**(-2)*0.5e-65,linestyle='--',color='k')
	#axs[1].annotate(r'$\propto k^{-2}$',(1e-17,1e-31))
	plt.xlim(k.min(),k.max())#2.7e-19,1e-16)
	axs[0].legend(frameon=False,fontsize=9)
	return fig,axs,k

def spectrum_graph(velfiles,boxsize):
	fig,axs=plt.subplots(2,sharex=True)
	marker=np.array(['v','.','+','2'])
	labels=np.array(['16 cells','32 cells','64 cells','128 cells'])
	for i in range(len(velfiles)):
		vx1,vy1,vz1,x,y,z=read_3cube(velfiles[i])
		vx1,vy1,vz1=vx1*code_units.v_cu, vy1*code_units.v_cu, vz1*code_units.v_cu
		x,y,z=x/x.max() * boxsize, y/y.max() * boxsize, z/z.max() * boxsize
		for j in range(2):
			if j==1:
				vx,vy,vz=subtract_radial(vx1,vy1,vz1,x,y,z,boxsize)
			if j==0:
				vx,vy,vz=vx1,vy1,vz1
			Ax=np.fft.fftn(vx)
			Ay=np.fft.fftn(vy)
			Az=np.fft.fftn(vz)
			A=np.sqrt(abs(Ax)**2+abs(Ay)**2+abs(Az)**2)
			Ncube=int(len(A[0,0,:]))
			k=np.fft.fftfreq(Ncube)*Ncube
			kx,ky,kz=np.meshgrid(k,k,k)
			K=np.sqrt(kx**2+ky**2+kz**2)
			bins=np.linspace(1,int(K.max()),int(K.max()))
			av1,ks1,args=binned_statistic(K.flatten(),abs(A/Ncube**3).flatten()**2,bins=bins)
			dk=ks1[1]-ks1[0]
			energy1=av1*4*np.pi*ks1[1:]**2 *dk
			if j==0:
				axs[j].loglog(ks1[:-1],energy1,marker[i],markersize=4,label=labels[i])
			if j==1:
				axs[j].loglog(ks1[:-1],energy1,marker[i],markersize=4)
			axs[j].set_ylabel(r'$P_{v}$',fontsize=11)	
			axs[j].tick_params(axis="y", labelsize=11,direction="in")
	axs[1].tick_params(axis="x", labelsize=11,direction="in")
	axs[1].set_xlabel(r'$cycles per box length$]',fontsize=11)
	plt.xlim(k.min(),k.max())
	axs[0].legend(frameon=False,fontsize=9)
	return fig,axs,k

	

