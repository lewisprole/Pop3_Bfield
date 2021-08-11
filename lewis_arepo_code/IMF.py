import numpy as np
import matplotlib.pyplot as plt
import arepo_utils 
import code_units 
import astropy.constants as ap
 

def energy_check(a,min_eject_length):
	masses=np.array([])
	ejects=0
	args=np.array([])
	x=sum(a.sinkmass*a.sinkx)/sum(a.sinkmass)
	y=sum(a.sinkmass*a.sinky)/sum(a.sinkmass)
	z=sum(a.sinkmass*a.sinkz)/sum(a.sinkmass)
	r_sinks=np.sqrt((a.sinkx-x)**2+(a.sinky-y)**2+(a.sinkz-z)**2) *code_units.d_cu
	r_cells=np.sqrt((a.x-x)**2+(a.y-y)**2+(a.z-z)**2) *code_units.d_cu
	for i in range(len(a.sinkmass)):
		if r_sinks[i]>(min_eject_length*code_units.d_cu):
			v=np.sqrt(a.sinkvx[i]**2+a.sinkvy[i]**2+a.sinkvz[i]**2) *code_units.v_cu
			mask=np.where(r_cells<r_sinks[i])
			mask_sinks=np.where(r_sinks<r_sinks[i])
			Mtot=sum(a.mass[mask])+sum(a.sinkmass[mask_sinks])

			Ek=(0.5*a.sinkmass[i]*code_units.M_cu*(v)**2)
			Eu= (ap.G.cgs.value * a.sinkmass[i]*code_units.M_cu * Mtot*code_units.M_cu) / r_sinks[i]
			if Ek>2*Eu:
				ejects+=1
				masses=np.append(masses,a.sinkmass[i])
				args=np.append(args,i)
	return args.astype(int)

def compiler(files):
	mass=np.array([])
	eject=np.array([])
	for i in range(len(files)):
		a=arepo_utils.aread(files[i])
		args=energy_check(a,0.03)
		mass=np.append(mass,a.sinkmass)
		eject=np.append(eject,a.sinkmass[args])
	return mass,eject

def IMF_time():
	fig,ax=plt.subplots(5,sharex=True)
	plt.subplots_adjust(hspace=0)

	colors='b','g','r','cyan','Purple'
	rhos=r'$\rho_{sink}$=10$^{-10}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-9}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-8}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-7}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-6}$gcm$^{-3}$'

	dirname='/scratch/c.c1521474/resolution_test/'

	A8=dirname+'/merge/1e8_projections/snapshot_135'
	A9=dirname+'/merge/1e9_projections/snapshot_135'
	A10=dirname+'/merge/1e10_projections/snapshot_135'
	A11=dirname+'/merge/1e11_projections/snapshot_132'
	A12=dirname+'/merge/1e12_projections/snapshot_495'

	B8=dirname+'/seed4/1e8/snapshot_191'
	B9=dirname+'/seed4/1e9/snapshot_191'
	B10=dirname+'/seed4/1e10/snapshot_191'
	B11=dirname+'/seed4/1e11/snapshot_190'
	B12=dirname+'/seed4/1e12/snapshot_488'

	C8=dirname+'/seed5/1e8/snapshot_166'
	C9=dirname+'/seed5/1e9/snapshot_166'
	C10=dirname+'/seed5/1e10/snapshot_166'
	C11=dirname+'/seed5/1e11/snapshot_166'
	C12=dirname+'/seed5/1e12/snapshot_445'

	mass8,eject8=compiler((A8,B8,C8))
	mass9,eject9=compiler((A9,B9,C9))
	mass10,eject10=compiler((A10,B10,C10))
	mass11,eject11=compiler((A11,B11,C11))
	mass12,eject12=compiler((A12,B12,C12))
	
	minmass=mass12.min()*0.7
	maxmass=mass8.max()*1.3
	bins=10**np.linspace(np.log10(minmass),np.log10(maxmass),50)
	mass=(mass8,mass9,mass10,mass11,mass12)
	eject=(eject8,eject9,eject10,eject11,eject12)
	for i in range(5):
		N,M=np.histogram(mass[i],bins=bins)
		ax[i].hist(mass[i],bins=bins,color=colors[i])
		ax[i].hist(eject[i],bins=bins,color='k')
		ax[i].set_ylim(0,N.max()+1)
		ax[i].text(1.22,0.5,rhos[i],ha='center', va='center', transform=ax[i].transAxes,fontsize=10)
	plt.subplots_adjust(left = 0.15,bottom = 0.17,right=0.7)
	ax[-1].set_xscale('log')
	ax[-1].set_xlabel(r'M [M$_\odot$yr$^{-1}$',fontsize=10)
	ax[2].set_ylabel(r'N$_{\rm sink}        $',fontsize=10,rotation=0)
