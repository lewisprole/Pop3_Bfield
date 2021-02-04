import arepo_utils 
from matplotlib.patches import Circle
from matplotlib.lines import Line2D
import numpy as np
import code_units
import astropy.constants as ap
from scipy.interpolate import interp1d
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic  

'''|||||||| simple plot showing T vs rho relationship |||||||||||||'''
def baro(filename):
	a=arepo_utils.aread(filename)
	x,y,z=np.histogram2d( np.log10(a.temp),np.log10(a.rho*code_units.rho_cu),bins=(800,800))
	plt.imshow(x/x,cmap='summer',aspect='auto',extent=[z[0],z[-1],y[-1],y[0]])
	plt.ylim(y[0],y[-1])
	plt.ylabel(r'log$_{10}$(T [k])',fontsize=10)
	plt.xlabel(r'log$_{10}$($\rho$ [g cm$^{-3}$])',fontsize=10)
	
	plt.tick_params(axis="x", labelsize=10,direction="in")
	plt.tick_params(axis="y", labelsize=10,direction="in")
	plt.subplots_adjust(left = 0.15,bottom = 0.17,right=0.9)





'''||||||functions to find volumetric heating rates||||||||||||'''

def compression(a):
	return ap.k_B.cgs.value * a.temp / ap.m_p.cgs.value * np.sqrt(32*ap.G.cgs.value/(3*np.pi)) * (a.rho*code_units.rho_cu)**(3/2)

def de_dt(a):
	heat=np.array([8,9,10,11,12,20,21,22,23,24,25])
	cool=np.array([0,1,2,3,4,6,7,13,14,15,16,17,18,19,26])
	heating = -np.sum(a.cooling[:,heat],1)
	cooling =np.sum(a.cooling[:,cool],1)
	acc=a.cooling[:,-1]
	return heating,cooling,acc

def cool_plot(file,AX):
	a=arepo_utils.aread(file)
	heating,cooling,acc =de_dt(a)
	comp=compression(a)


	#COOLING,y,z=np.histogram2d( np.log10(cooling)[~np.isnan(np.log10(cooling))],  np.log10(a.rho*code_units.rho_cu)[~np.isnan(np.log10(cooling))],  bins=(500,500))
	#AX.imshow(COOLING/COOLING,cmap='winter',aspect='auto',label=r'$\Lambda$',extent=[z[0],z[-1],y[-1],y[0]])
	COOLING,rho,z=binned_statistic(a.rho*code_units.rho_cu,cooling,bins = 10**np.linspace(np.log10(a.rho.min()*code_units.rho_cu),np.log10(a.rho.max()*code_units.rho_cu),50))
	AX.loglog(rho[:-1],COOLING,'cyan',linewidth=5)

	#HEATING,y,z=np.histogram2d( np.log10(heating)[~np.isnan(np.log10(heating))],np.log10(a.rho*code_units.rho_cu)[~np.isnan(np.log10(heating))],bins=(500,500))
	#AX.imshow(HEATING/HEATING,cmap='RdYlBu',aspect='auto',label=r'$-\Gamma$',extent=[z[0],z[-1],y[-1],y[0]])
	HEATING,rho,z=binned_statistic(a.rho*code_units.rho_cu,heating,bins = 10**np.linspace(np.log10(a.rho.min()*code_units.rho_cu),np.log10(a.rho.max()*code_units.rho_cu),50))
	AX.loglog(rho[:-1],HEATING,'r')

	#COMP,y,z=np.histogram2d( np.log10(comp),np.log10(a.rho*code_units.rho_cu),bins=(500,500))
 	#AX.imshow(COMP/COMP,cmap='summer',aspect='auto',label=r'$-\Gamma$',extent=[z[0],z[-1],y[-1],y[0]])
	COMP,rho,z=binned_statistic(a.rho*code_units.rho_cu,comp,bins = 10**np.linspace(np.log10(a.rho.min()*code_units.rho_cu),np.log10(a.rho.max()*code_units.rho_cu),50))
	AX.loglog(rho[:-1],COMP,'k')

	#ACC,y,z=np.histogram2d( np.log10(acc)[~np.isnan(np.log10(acc))],np.log10(a.rho*code_units.rho_cu)[~np.isnan(np.log10(acc))],bins=(500,500))
	#AX.imshow(ACC/ACC,cmap='spring',aspect='auto',label=r'$Gamma_L$',extent=[z[0],z[-1],y[-1],y[0]])
	ACC,rho,z=binned_statistic(a.rho*code_units.rho_cu,acc,bins = 10**np.linspace(np.log10(a.rho.min()*code_units.rho_cu),np.log10(a.rho.max()*code_units.rho_cu),50))
	AX.loglog(rho[:-1],ACC,'fuchsia')
	
	#AX.set_ylim(y[0],y[-1])
	AX.tick_params(axis="x", labelsize=9,direction="in")
	AX.tick_params(axis="y", labelsize=9,direction="in")
	#AX.set_ylim(y[0],y[-1])
	#return z[0],z[-1]

def pannel_plot(file8,file9,file10,file11,file12):
	labels=r'$\Gamma$',r'$\Gamma_L$',r'$\frac{k_B T}{m_p} \sqrt{\frac{32G}{3\pi}}$($\rho$)$^{3/2}$',r'$-\Lambda$'
	fig,axs=plt.subplots(5,sharex=True)
	plt.subplots_adjust(hspace=0,top=0.95,bottom=0.12,right=0.75,left=0.15)

	#X1,X2=
	cool_plot(file12,axs[4])
	axs[4].set_ylim(10**-19,10**19)
	
	axs[4].text(0.18,0.9,r'$\rho_{\rm sink}$=10$^{-6}$gcm$^{-3}$',ha='center', va='center', transform=axs[4].transAxes,fontsize=10)

	#x1,x2=
	cool_plot(file11,axs[3])
	axs[3].set_ylim(10**-19,10**19)
	axs[3].text(0.18,0.88,r'$\rho_{\rm sink}$=10$^{-7}$gcm$^{-3}$',ha='center', va='center', transform=axs[3].transAxes,fontsize=10)

	#x1,x2=
	cool_plot(file10,axs[2])
	axs[2].set_ylim(10**-19,10**19)
	axs[2].text(0.18,0.88,r'$\rho_{\rm sink}$=10$^{-8}$gcm$^{-3}$',ha='center', va='center', transform=axs[2].transAxes,fontsize=10)

	#x1,x2=
	cool_plot(file9,axs[1])
	axs[1].set_ylim(10**-19,10**19)
	axs[1].text(0.18,0.88,r'$\rho_{\rm sink}$=10$^{-9}$gcm$^{-3}$',ha='center', va='center', transform=axs[1].transAxes,fontsize=10)

	#x1,x2=
	cool_plot(file8,axs[0])
	axs[0].set_ylim(10**-19,10**19)
	axs[0].text(0.18,0.88,r'$\rho_{\rm sink}$=10$^{-10}$gcm$^{-3}$',ha='center', va='center', transform=axs[0].transAxes,fontsize=10)
	axs[0].set_xlim(10**-16,10**-3.5)


	axs[0].axvline(x=1e8*rho_cu,linestyle='--',color='k')
	axs[1].axvline(x=1e9*rho_cu,linestyle='--',color='k')
	axs[2].axvline(x=1e10*rho_cu,linestyle='--',color='k')
	axs[3].axvline(x=1e11*rho_cu,linestyle='--',color='k')
	axs[4].axvline(x=1e12*rho_cu,linestyle='--',color='k')

	#plt.subplots_adjust(top=0.95,bottom=0.05)
	axs[4].set_xlabel(r'Log$_{10}(\rho$ [gcm$^{-3}$])',fontsize=10)
	axs[2].set_ylabel(r'Log$_{10}$($\frac{\rm de}{\rm dt}$ [erg s$^{-1}$ cm$^{-3}]$)')
	axs[4].set_xlim(10**(-18),10**(-3.5))
	line1=Line2D([0], [0], color='r')#linestyle='none', marker='o', markerfacecolor='r',markeredgecolor='r')
	line2=Line2D([0], [0], color='fuchsia')#linestyle='none', marker='o', markerfacecolor='fuchsia',markeredgecolor='fuchsia')
	line3=Line2D([0], [0], color='k')#linestyle='none', marker='o', markerfacecolor='k',markeredgecolor='k')
	line4=Line2D([0], [0], color='cyan')#linestyle='none', marker='o',markerfacecolor='cyan',markeredgecolor='cyan')
	axs[0].legend([line1,line2,line3,line4],(r'$\Gamma$',r'$\Gamma_L$',r'$\frac{k_B T}{m_p} \sqrt{\frac{32G}{3\pi}}$($\rho$)$^{3/2}$',r'$-\Lambda$'),fontsize=10,frameon=False,markerscale=1,loc=(0.99,0.1))


'''|||||||||| functions to plot abundances vs density||||||||||'''
def abundance(files):
	colors='r','k','cyan','pink','green','b'
	labels=r'H$_{\rm 2}$',r'H$^{+}$',r'D$^{+}$','HD',r'He$^{+}$',r'He$^{++}$'
	texts=r'$\rho_{\rm sink}$=10$^{-10}$gcm$^{-3}$',r'$\rho_{\rm sink}$=10$^{-9}$gcm$^{-3}$',r'$\rho_{\rm sink}$=10$^{-8}$gcm$^{-3}$',r'$\rho_{\rm sink}$=10$^{-7}$gcm$^{-3}$',r'$\rho_{\rm sink}$=10$^{-6}$gcm$^{-3}$'
	fig,axs=plt.subplots(5,sharex=True)
	plt.subplots_adjust(hspace=0)
	for i in range(len(files)):
		a=arepo_utils.aread(files[i])
		for j in range(a.chem.shape[1]):
			abund,rho,z=binned_statistic(a.rho*code_units.rho_cu,a.chem[:,j],bins = 10**np.linspace(np.log10(a.rho.min()*code_units.rho_cu),np.log10(a.rho.max()*code_units.rho_cu),50))
			if i==0:
				axs[i].loglog(rho[:-1],abund,colors[j],label=labels[j])
			else:
				axs[i].loglog(rho[:-1],abund,colors[j])
			
			axs[i].set_ylim(1e-49,1e10)
 			axs[i].tick_params(axis="x", labelsize=9,direction="in")
			axs[i].tick_params(axis="y", labelsize=9,direction="in")

		axs[i].text(0.2,0.1,texts[i],ha='center', va='center', transform=axs[i].transAxes,fontsize=10)
		axs[i].set_yticks([1e-41,1e-26,1e-11,1e4])

	axs[4].set_xlim(1e-16,10**-3.5)
	axs[0].axvline(x=1e8*code_units.rho_cu,linestyle='--',color='k')
	axs[1].axvline(x=1e9*code_units.rho_cu,linestyle='--',color='k')
	axs[2].axvline(x=1e10*code_units.rho_cu,linestyle='--',color='k')
	axs[3].axvline(x=1e11*code_units.rho_cu,linestyle='--',color='k')
	axs[4].axvline(x=1e12*code_units.rho_cu,linestyle='--',color='k')
	axs[4].set_xlabel(r'Log$_{10}(\rho$ [gcm$^{-3}$])',fontsize=10)
	axs[2].set_ylabel(r'$X$',fontsize=10)
	plt.subplots_adjust(hspace=0,top=0.95,bottom=0.12,right=0.75,left=0.15)
	axs[0].legend(fontsize=10,frameon=False,loc=(1.01,-0.25))


'''||||||FUNCTIONS NOT USED||||||||'''



'''||||||| average mass to calculate free-fall time  |||||||'''


def Tff(a):
	'''calculates cumulative mass within shells and used it to calculate free-fall time'''
	mid=np.where(a.rho==a.rho.max())
	r=np.sqrt((a.x-a.x[mid])**2+(a.y-a.y[mid])**2+(a.z-a.z[mid])**2) *code_units.d_cu
	args=np.argsort(r)
	rs=10**np.linspace(np.log10(np.sort(r)[1]),np.log10(r.max()),100)
	M=np.zeros_like(rs)
	for i in range(len(rs)-1):
		M[i]=np.sum(a.mass[np.where(r<rs[i+1])])*code_units.M_cu
	M[-1]=M[-2]
	rho=M/(4/3 *np.pi * rs**3)
	tff=np.sqrt(3/(32*ap.G.cgs.value*rho))
	return tff,rs

def function(tff,rs):
	'''interpolate to fit free-fall time back to original data'''
	f=interp1d(rs,tff,fill_value="extrapolate")
	return f




