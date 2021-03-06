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
 	HEATING,y,z=np.histogram2d( np.log10(heating)[~np.isnan(np.log10(heating))],np.log10(a.rho*code_units.rho_cu)[~np.isnan(np.log10(heating))],bins=(500,500))
 	AX.imshow(HEATING/HEATING,cmap='RdYlBu',aspect='auto',label=r'$-\Gamma$',extent=[z[0],z[-1],y[-1],y[0]])

	COMP,y,z=np.histogram2d( np.log10(comp),np.log10(a.rho*code_units.rho_cu),bins=(500,500))
 	AX.imshow(COMP/COMP,cmap='summer',aspect='auto',label=r'$-\Gamma$',extent=[z[0],z[-1],y[-1],y[0]])

	ACC,y,z=np.histogram2d( np.log10(acc)[~np.isnan(np.log10(acc))],np.log10(a.rho*code_units.rho_cu)[~np.isnan(np.log10(acc))],bins=(500,500))
	AX.imshow(ACC/ACC,cmap='spring',aspect='auto',label=r'$Gamma_L$',extent=[z[0],z[-1],y[-1],y[0]])

	COOLING,y,z=np.histogram2d( np.log10(cooling)[~np.isnan(np.log10(cooling))],  np.log10(a.rho*code_units.rho_cu)[~np.isnan(np.log10(cooling))],  bins=(500,500))
	AX.imshow(COOLING/COOLING,cmap='winter',aspect='auto',label=r'$\Lambda$',extent=[z[0],z[-1],y[-1],y[0]])
	
	AX.set_ylim(y[0],y[-1])
	AX.tick_params(axis="x", labelsize=9,direction="in")
	AX.tick_params(axis="y", labelsize=9,direction="in")
	AX.set_ylim(y[0],y[-1])
	return z[0],z[-1]

def pannel_plot(file9,file10,file11,file12):
	fig,axs=plt.subplots(5,sharex=True)
	plt.subplots_adjust(hspace=0,top=0.95,bottom=0.1,right=0.8,left=0.15)

	X1,X2=cool_plot(file12,axs[4])
	axs[4].set_ylim(-16,13)
	axs[4].text(0.82,0.1,r'$\rho_{\rm sink}$=10$^{-6}$gcm$^{-3}$',ha='center', va='center', transform=axs[4].transAxes,fontsize=10)

	x1,x2=cool_plot(file11,axs[3])
	axs[3].set_ylim(-16,13)
	axs[3].text(0.82,0.1,r'$\rho_{\rm sink}$=10$^{-7}$gcm$^{-3}$',ha='center', va='center', transform=axs[3].transAxes,fontsize=10)

	x1,x2=cool_plot(file10,axs[2])
	axs[2].set_ylim(-16,13)
	axs[2].text(0.82,0.1,r'$\rho_{\rm sink}$=10$^{-8}$gcm$^{-3}$',ha='center', va='center', transform=axs[2].transAxes,fontsize=10)

	x1,x2=cool_plot(file9,axs[1])
	axs[1].set_ylim(-16,13)
	axs[1].text(0.82,0.1,r'$\rho_{\rm sink}$=10$^{-9}$gcm$^{-3}$',ha='center', va='center', transform=axs[1].transAxes,fontsize=10)

	x1,x2=cool_plot(file8,axs[0])
	axs[0].set_ylim(-16,13)
	axs[0].text(0.82,0.1,r'$\rho_{\rm sink}$=10$^{-10}$gcm$^{-3}$',ha='center', va='center', transform=axs[0].transAxes,fontsize=10)

	plt.subplots_adjust(top=0.95,bottom=0.05)
	axs[4].set_xlabel(r'Log$_{10}(\rho$ [gcm$^{-3}$])',fontsize=10)
	axs[2].set_ylabel(r'Log$_{10}$($\frac{\rm de}{\rm dt}$ [erg s$^{-1}$ cm$^{-3}]$)')
	axs[4].set_xlim(-18,-4)
	line1=Line2D([0], [0], linestyle='none', marker='o', markerfacecolor='brown',markeredgecolor='brown')
	line2=Line2D([0], [0], linestyle='none', marker='o', markerfacecolor='fuchsia',markeredgecolor='fuchsia')
	line3=Line2D([0], [0], linestyle='none', marker='o', markerfacecolor='forestgreen',markeredgecolor='forestgreen')
	line4=Line2D([0], [0], linestyle='none', marker='o',markerfacecolor='b',markeredgecolor='b')
	axs[0].legend([line1,line2,line3,line4],(r'$\Gamma$',r'$\Gamma_L$',r'$\frac{k_B T}{m_p} \sqrt{\frac{32G}{3\pi}}$($\rho$)$^{3/2}$',r'$-\Lambda$'),fontsize=10,frameon=False,markerscale=1,loc=(0.99,0.1))









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








''' ||||||||||| compressional heating rate per mass |||||||||||'''
def compression(a):
	return ap.k_B.cgs.value * a.temp / ap.m_p.cgs.value * np.sqrt(32*ap.G.cgs.value/(2*np.pi)) * np.sqrt(a.rho*code_units.rho_cu)
def compression_line(a):
	rhomin,rhomax=np.log10(a.rho.min()*code_units.rho_cu),np.log10(a.rho.max()*code_units.rho_cu)
	comp,rhos,z=binned_statistic(a.rho*code_units.rho_cu,ap.k_B.cgs.value * a.temp / ap.m_p.cgs.value * np.sqrt(32*ap.G.cgs.value/(2*np.pi)) * np.sqrt(a.rho*code_units.rho_cu),bins=10**np.linspace(rhomin,rhomax,100))
	return comp,rhos[1:]


'''||||||||| working out the accretion heating rate per mass |||||||||||'''

def accretion(file1,file2,Rstar):
	a=arepo_utils.aread(file2)
	b=arepo_utils.aread(file1)
	dt=(a.time-b.time )*code_units.t_cu
	heating = np.zeros_like(a.rho)
	for i in range(len(a.sinkmass)):
		dM=(a.sinkmass[i]-b.sinkmass[i])*code_units.M_cu
		rate=dM/dt
		L=ap.G.cgs.value*a.sinkmass[i]*code_units.M_cu * rate / Rstar
		r=np.sqrt((a.x-a.sinkx[i])**2+(a.y-a.sinky[i])**2+(a.z-a.sinkz[i])**2) *code_units.d_cu
		heating +=  a.rho*code_units.rho_cu * L /(4*np.pi * r**2)   #per mass 
	return heating 
		

''' |||||||| e = rho * u, de/dt = (drho/dt)u + (du/dt)rho, we want du/dt (rate per mass)'''
'''---> du/dt = [de/dt - (drho/dt)u]/rho, use continuity equation for drho/dt |||||||'''
	
def du_dt(a):
	drho_dt =  -a.divv*(1/code_units.t_cu) * a.rho*code_units.rho_cu
	
	#create args that separate heating and cooling terms in a.cooling
	heat=np.array([8,9,10,11,12,20,21,22,23,24,25])
	cool=np.array([0,1,2,3,4,6,7,13,14,15,16,17,18,19,26])
	
	#convert cooling rate per volume into cooling rate per mass 
	for i in range(a.cooling.shape[1]):
		a.cooling[:,i] =  (a.cooling[:,i] - drho_dt* a.u*code_units.v_cu**2) / (a.rho*code_units.rho_cu)  #a.cooling[:,i] - drho_dt

	#total heating and  cooling rate per mass
	heat_rate= np.sum(a.cooling[:,heat],1) 
	cool_rate= np.sum(a.cooling[:,cool],1) 
	
	return heat_rate, cool_rate 

'''|||||||plot heating and cooling vs density|||||||'''


def plot(file,Rstar):
	a=arepo_utils.aread(file)
	heating,cooling =du_dt(a)
	comp=compression_line(a)

	HEATING,y,z=np.histogram2d( np.log10(-heating)[~np.isnan(np.log10(-heating))],       np.log10(a.rho*code_units.rho_cu)[~np.isnan(np.log10(-heating))],bins=(800,800))
	plt.imshow(HEATING/HEATING,cmap='RdYlBu',aspect='auto',label=r'$-\Gamma$',extent=[z[0],z[-1],y[-1],y[0]])

	COMP,y,z=np.histogram2d( np.log10(comp),np.log10(a.rho*code_units.rho_cu),bins=(800,800))
	plt.imshow(COMP/COMP,cmap='spring',aspect='auto',label=r'$k_B T/m_p \sqrt{\frac{32G\rho}{3\pi}$',extent=[z[0],z[-1],y[-1],y[0]])

	COOLING,y,z=np.histogram2d( np.log10(cooling)[~np.isnan(np.log10(cooling))],  np.log10(a.rho*code_units.rho_cu)[~np.isnan(np.log10(cooling))],  bins=(800,800))
	plt.imshow(COOLING/COOLING,cmap='winter',aspect='auto',label=r'$\Lambda$',extent=[z[0],z[-1],y[-1],y[0]])

        #TOTHEAT,y,z=np.histogram2d( np.log10(-heating+comp+acc),np.log10(a.rho*code_units.rho_cu),bins=(800,800))
        #plt.imshow(TOTHEAT/TOTHEAT,cmap='afmhot',aspect='auto',label=r'$\Lambda$',extent=[z[0],z[-1],y[-1],y[0]])

	plt.ylim(y[0],y[-1])

	plt.ylabel(r'erg s$^{-1}$ g $^{-1}$')
	plt.xlabel(r'$\rho$ [g cm$^{-3}$]')

	line1=Line2D([0], [0], color='brown', lw=2)
	line3=Line2D([0], [0], color='forestgreen', lw=2)
	line4=Line2D([0], [0], color='b', lw=2)

	plt.legend([line1,line3,line4],(r'$\Gamma$',r'$\frac{k_B T}{m_p} \sqrt{\frac{32G\rho}{3\pi}}$',r'$-\Lambda$'),fontsize=10,frameon=False,markerscale=10,loc='upper left')

def plot_sinks(file1,file2,Rstar):
	a=arepo_utils.aread(file2)
	heating,cooling =du_dt(a)
	comp=compression(a)
	acc=accretion(file1,file2,Rstar)

	HEATING,y,z=np.histogram2d( np.log10(-heating)[~np.isnan(np.log10(-heating))],       np.log10(a.rho*code_units.rho_cu)[~np.isnan(np.log10(-heating))],bins=(800,800))
	plt.imshow(HEATING/HEATING,cmap='RdYlBu',aspect='auto',label=r'$-\Gamma$',extent=[z[0],z[-1],y[-1],y[0]])

	ACC,y,z=np.histogram2d( np.log10(acc),np.log10(a.rho*code_units.rho_cu),bins=(800,800))
	plt.imshow(ACC/ACC,cmap='spring',aspect='auto',label=r'$Gamma_L$',extent=[z[0],z[-1],y[-1],y[0]])

	COMP,y,z=np.histogram2d( np.log10(comp),np.log10(a.rho*code_units.rho_cu),bins=(800,800))
	plt.imshow(COMP/COMP,cmap='summer',aspect='auto',label=r'$k_B T/m_p \sqrt{\frac{32G\rho}{3\pi}$',extent=[z[0],z[-1],y[-1],y[0]])

	COOLING,y,z=np.histogram2d( np.log10(cooling)[~np.isnan(np.log10(cooling))],  np.log10(a.rho*code_units.rho_cu)[~np.isnan(np.log10(cooling))],  bins=(800,800))
	plt.imshow(COOLING/COOLING,cmap='winter',aspect='auto',label=r'$\Lambda$',extent=[z[0],z[-1],y[-1],y[0]])

	#TOTHEAT,y,z=np.histogram2d( np.log10(-heating+comp+acc),np.log10(a.rho*code_units.rho_cu),bins=(800,800))
	#plt.imshow(TOTHEAT/TOTHEAT,cmap='afmhot',aspect='auto',label=r'$\Lambda$',extent=[z[0],z[-1],y[-1],y[0]])

	plt.ylim(y[0],y[-1])

	plt.ylabel(r'log$_{10}$(erg s$^{-1}$ g $^{-1}$)',fontsize=10)
	plt.xlabel(r'log$_{10}$($\rho$ [g cm$^{-3}$])',fontsize=10)

	line1=Line2D([0], [0], color='brown', lw=2)
	line2=Line2D([0], [0], color='fuchsia', lw=2)
	line3=Line2D([0], [0], color='forestgreen', lw=2)
	line4=Line2D([0], [0], color='b', lw=2)
	plt.legend([line1,line2,line3,line4],(r'$\Gamma$',r'$\Gamma_L$',r'$\frac{k_B T}{m_p} \sqrt{\frac{32G\rho}{3\pi}}$',r'$-\Lambda$'),fontsize=10,frameon=False,markerscale=10,loc='upper left')
	plt.tick_params(axis="x", labelsize=10,direction="in")
	plt.tick_params(axis="y", labelsize=10,direction="in")
	plt.subplots_adjust(left = 0.15,bottom = 0.17,right=0.9)


def prep_panel(file1,file2,Rstar,AX):
	a=arepo_utils.aread(file2)
	heating,cooling =du_dt(a)
	comp,rhos=compression_line(a)
	acc=accretion(file1,file2,Rstar)
	HEATING,y,z=np.histogram2d( np.log10(-heating)[~np.isnan(np.log10(-heating))],       np.log10(a.rho*code_units.rho_cu)[~np.isnan(np.log10(-heating))],bins=(100,500))
	AX.imshow(HEATING/HEATING,cmap='RdYlBu',aspect='auto',label=r'$-\Gamma$',extent=[z[0],z[-1],y[-1],y[0]])
	AX.plot(np.log10(rhos),np.log10(comp),c='indigo')
	ACC,y,z=np.histogram2d( np.log10(acc),np.log10(a.rho*code_units.rho_cu),bins=(100,500))
	AX.imshow(ACC/ACC,cmap='spring',aspect='auto',label=r'$Gamma_L$',extent=[z[0],z[-1],y[-1],y[0]])
	COOLING,y,z=np.histogram2d( np.log10(cooling)[~np.isnan(np.log10(cooling))],  np.log10(a.rho*code_units.rho_cu)[~np.isnan(np.log10(cooling))],  bins=(100,500))
	AX.imshow(COOLING/COOLING,cmap='winter',aspect='auto',label=r'$\Lambda$',extent=[z[0],z[-1],y[-1],y[0]])
	AX.set_ylim(y[0],y[-1])
	AX.tick_params(axis="x", labelsize=9,direction="in")
	AX.tick_params(axis="y", labelsize=9,direction="in")
	return z[0],z[-1]

def plot_panels(file1_9,file2_9,file1_11,file2_11,file1_12,file2_12,Rstars):
	fig,axs=plt.subplots(5,sharex=True)	
	plt.subplots_adjust(hspace=0,top=0.95,bottom=0.1,right=0.8,left=0.15)
	X1,X2=prep_panel(file1_12,file2_12,Rstars[4],axs[4])
	axs[4].set_xlim(X1,X2)
	axs[4].set_ylim(-16,13)
	axs[4].text(0.82,0.1,r'$\rho_{\rm sink}$=10$^{-6}$gcm$^{-3}$',ha='center', va='center', transform=axs[4].transAxes,fontsize=10)
	x1,x2=prep_panel(file1_11,file2_11,Rstars[3],axs[3])
	axs[3].set_xlim(X1,X2)
	axs[3].set_ylim(-16,13)
	axs[3].text(0.82,0.1,r'$\rho_{\rm sink}$=10$^{-7}$gcm$^{-3}$',ha='center', va='center', transform=axs[3].transAxes,fontsize=10)
	x1,x2=prep_panel(file1_10,file2_10,Rstars[2],axs[2])
	axs[2].set_xlim(X1,X2)
	axs[2].set_ylim(-16,13)
	axs[2].text(0.82,0.1,r'$\rho_{\rm sink}$=10$^{-8}$gcm$^{-3}$',ha='center', va='center', transform=axs[2].transAxes,fontsize=10)
	x1,x2=prep_panel(file1_9,file2_9,Rstars[1],axs[1])
	axs[1].set_xlim(X1,X2)
	axs[1].set_ylim(-16,13)
	axs[1].text(0.82,0.1,r'$\rho_{\rm sink}$=10$^{-9}$gcm$^{-3}$',ha='center', va='center', transform=axs[1].transAxes,fontsize=10)
	x1,x2=prep_panel(file1_8,file2_8,Rstars[0],axs[0])
	axs[0].set_xlim(X1,X2)
	axs[0].set_ylim(-16,13)
	axs[0].text(0.82,0.1,r'$\rho_{\rm sink}$=10$^{-10}$gcm$^{-3}$',ha='center', va='center', transform=axs[0].transAxes,fontsize=10)
	plt.subplots_adjust(top=0.95,bottom=0.05)
	axs[4].set_xlabel(r'Log$_{10}(\rho$ [gcm$^{-3}$])',fontsize=10)
	axs[2].set_ylabel(r'Log$_{10}$($\frac{\rm du}{\rm dt}$ [erg s$^{-1}$ g$^{-1}]$)')
	axs[4].set_xlim(-18,-4)
	line1=Line2D([0], [0], linestyle='none', marker='o', markerfacecolor='brown',markeredgecolor='brown')
	line2=Line2D([0], [0], linestyle='none', marker='o', markerfacecolor='fuchsia',markeredgecolor='fuchsia')
	line3=Line2D([0], [0], linestyle='none', marker='o', markerfacecolor='forestgreen',markeredgecolor='forestgreen')
	line4=Line2D([0], [0], linestyle='none', marker='o',markerfacecolor='b',markeredgecolor='b')
	axs[0].legend([line1,line2,line3,line4],(r'$\Gamma$',r'$\Gamma_L$',r'$\frac{k_B T}{m_p} \sqrt{\frac{32G\rho}{3\pi}}$',r'$-\Lambda$'),fontsize=10,frameon=False,markerscale=1,loc=(0.99,0.1))


'''||| plot heating and cooling as functions of density and radial distance |||'''


def both(file):
	a=arepo_utils.aread(file)
	fig,axs=plt.subplots(nrows=1, ncols=2,sharey=True)
        #radial distances fro CoM
	mid=np.where(a.rho==a.rho.max())
	r=np.sqrt((a.x-a.x[mid])**2+(a.y-a.y[mid])**2+(a.z-a.z[mid])**2)

        #create args that separate heating and cooling terms in a.cooling
	heat=np.array([8,9,10,11,12,20,21,22,23,24,25])
	cool=np.array([0,1,2,3,4,6,7,13,14,15,16,17,18,19,26])

	#PLOT 1 

        #calculate free-fall times
	tff=np.sqrt(3*np.pi/(32*ap.G.cgs.value*a.rho*code_units.rho_cu))

        #2d histogram of heating vs rho
	t_heat= abs(a.u*code_units.v_cu**2 *a.rho*code_units.rho_cu / np.sum(a.cooling[:,heat],1) /(60*60*24*365))
	t_cool= abs(a.u*code_units.v_cu**2 *a.rho*code_units.rho_cu / np.sum(a.cooling[:,cool],1) /(60*60*24*365))

	x,y,z=np.histogram2d( np.log10(t_heat),np.log10(a.rho*code_units.rho_cu),bins=(800,800))
	line1=axs[0].imshow(x/x,cmap='RdYlBu',aspect='auto',label=r'$t_h$',extent=[z[0],z[-1],y[-1],y[0]])
	x,y,z=np.histogram2d( np.log10(t_cool),np.log10(a.rho*code_units.rho_cu),bins=(800,800))
	line2=axs[0].imshow(x/x,cmap='winter',aspect='auto',label=r'$t_c$',extent=[z[0],z[-1],y[-1],y[0]])
	axs[0].set_ylim(y[0],y[-1])

	rho=10**np.linspace(z[0],z[-1],100)
	tff=np.sqrt(3*np.pi/(32*ap.G.cgs.value*rho))/(60*60*24*365)
	#compression_rate =ap.k_B.cgs.value*a.temp/(ap.m_p.cgs.value) * np.sqrt(32*ap.G.cgs.value/(3*np.pi)) * (a.rho*code_units.rho_cu)**(3/2)
	#t_comp=(a.u*code_units.v_cu**2 *a.rho*code_units.rho_cu / compression_rate) /(60*60*24*365)
	line3=axs[0].plot(np.log10(rho),np.log10(tff),'k',label=r'$t_{ff}$')
	#line3=axs[0].plot(np.log10(a.rho*code_units.rho_cu),np.log10(t_comp),'k',label=r'$t_{ff}$')	

        #create legend patches
	line1=Line2D([0], [0], color='brown', lw=2)
	line2=Line2D([0], [0], color='b', lw=2)
	line3=Line2D([0], [0], color='k', lw=2)
	#axs[0].legend([line1,line2,line3],(r'$t_h$',r'$t_c$',r'$t_{ff}$'),fontsize=12,frameon=False,markerscale=10)

        #configure subplots
        #plt.subplots_adjust(left = 0.2,bottom = 0.17,right=0.9)
	axs[0].tick_params(axis="x", labelsize=9,direction="in")
	axs[0].tick_params(axis="y", labelsize=9,direction="in")
	axs[0].set_xlabel(r'log$_{10}(\rho$ [gcm$^{-3}$])',fontsize=10)
	axs[0].set_ylabel(r'log$_{10}$(t [yrs])',fontsize=10,rotation=90)
	axs[0].set_ylim(np.log10(t_cool).min()-1,y[-1])
	axs[0].set_xticks([-18,-14,-10,-6,-22])


	#PLOT 2 
	tff,rs=Tff(a)
	mask=np.where(r>0)
	#2d histogram of cooling vs rho
	x,y,z=np.histogram2d(np.log10(t_heat[mask]),np.log10(r[mask]*code_units.d_cu),bins=(800,800))
	line4=axs[1].imshow(x/x,cmap='RdYlBu',aspect='auto',label=r'$t_h$',extent=[z[0],z[-1],y[-1],y[0]])
	x,y,z=np.histogram2d(np.log10(t_cool[mask]),np.log10(r[mask]*code_units.d_cu),bins=(800,800))
	line5=axs[1].imshow(x/x,cmap='winter',aspect='auto',label=r'$t_c$',extent=[z[0],z[-1],y[-1],y[0]])
	line6=axs[1].plot(np.log10(rs),np.log10(tff/(60*60*24*365)),'k',label=r'$t_{ff}$')
	line1=Line2D([0], [0], color='brown', lw=2)
	line2=Line2D([0], [0], color='b', lw=2)
	line3=Line2D([0], [0], color='k', lw=2)
	axs[1].legend([line1,line2,line3],(r'$t_h$',r'$t_c$',r'$t_{ff}$'),fontsize=10,frameon=False,markerscale=10,loc='upper left')
	plt.subplots_adjust(left = 0.15,bottom = 0.252,top=0.85,right=0.95,wspace = 0)
	axs[1].tick_params(axis="x", labelsize=9,direction="in")
	axs[1].tick_params(axis="y", labelsize=9,direction="in")
	axs[1].set_xlabel(r'log$_{10}$(R [cm])',fontsize=10)
	axs[1].set_xticks([13,15,17,19])
 	#axs[1].plt.ylabel(r'log$_{10}(t$ [yrs])       ',fontsize=20,rotation=90)
