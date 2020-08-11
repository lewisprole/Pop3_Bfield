import numpy as np
import matplotlib.pyplot as plt
import code_units
import sys
import io
import arepo_utils
import astropy.constants as ap


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
			mask=np.where(rs<x[i+1])
		else:
			mask=np.where((rs>x[i]) & (rs<x[i+1]))
		RS=rs[mask]
		if sum(weight[mask])>0:
			av=sum(variable[mask]*weight[mask])/sum(weight[mask])
		else:
			av=0
		avs=np.append(avs,av)
	return x,avs
	
def cycle_plot(snaps,weight_type,bins,labels):
	
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


