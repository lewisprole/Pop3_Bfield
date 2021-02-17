
import numpy as np
import matplotlib.pyplot as plt 
import arepo_utils
import io 
import sys
import code_units
import astropy.constants as ap
import os
from matplotlib.lines import Line2D
from scipy.stats import binned_statistic
from scipy.interpolate import interp1d


def snapname(start,i,interval):
	'''creates snapshot id'''
	n='00'+str(start+i*interval)
	if start+i*interval>9:
		n='0'+str(start+i*interval)
	if start+i*interval>99:
		n=str(start+i*interval)
	return n


def cycle(dirname,start,end,interval,name):
	if os.path.isfile(name):
		print("file exists, appending")
		f=open(name, "a+")
	else:
		print("creating new file")
		f = open(name, "x")
	Mtot=[]
	N=[]
	t=[]
	num=int((end-start)/interval)
	for i in range(num):
		text_trap = io.StringIO() #prevent massive text output from snapshot reads
		sys.stdout = text_trap
		n=snapname(start,i,interval)
		a=arepo_utils.aread(dirname+'snapshot_'+n)
		sys.stdout = sys.__stdout__
		if a.npart[-1]>0:
			Mtot.append(sum(a.sinkmass))
			N.append(a.npart[-1])
			t.append(a.time)
			f.write(str(a.npart[-1]) + ' ' + str(sum(a.sinkmass)) + ' ' + str(a.time) + '\n')
		print(str(n)+' :done')
	f.close()
	return Mtot,N,t
	

def txtread(txtfile):
	N=[]
	M=[]
	t=[]
	with open(txtfile) as f:
		for line in f.readlines():
			N.append(line.split()[0])
			M.append(line.split()[1])
			t.append(line.split()[2])
	return np.asarray(N).astype(float),np.asarray(M).astype(float) ,np.asarray(t).astype(float)

def plot_MN(files):
	Nmax=0
	fig,axs=plt.subplots(2,sharex=True)
	plt.subplots_adjust(wspace=0, hspace=0)
	colors='b','g','r','cyan','purple'
	labels=r'$\rho_{sink}$=10$^{-10}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-9}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-8}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-7}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-6}$gcm$^{-3}$'
	for i in range(len(files)):
		N,M,t=txtread(files[i])
		if i==0:
			t0=t[0]
		axs[0].plot((t-t0)*code_units.t_cu/(60*60*24*365),N,color=colors[i],label=labels[i])
		axs[1].plot((t-t0)*code_units.t_cu/(60*60*24*365),M*code_units.M_cu/ap.M_sun.cgs.value,color=colors[i],label=labels[i])
		if N.max()>Nmax:
			Nmax=N.max()

	axs[1].set_xlabel(r'$t \ [yrs]$',fontsize=10)
	axs[0].set_ylabel(r'$N_{sinks}$',fontsize=10)
	axs[1].set_ylabel(r'$\sum M_{sink} \ [M_{\odot}]$',fontsize=10)
	axs[0].set_yticks(np.arange(1,Nmax+2,2))
	axs[1].tick_params(axis="x", labelsize=10,direction="in")
	axs[0].tick_params(axis="x", labelsize=10,direction="in")
	axs[0].tick_params(axis="y", labelsize=10,direction="in")
	axs[1].tick_params(axis="y", labelsize=10,direction="in")
	axs[1].legend(fontsize=10,loc='lower right',frameon=False)#,bbox_to_anchor=(0.99, 1))
	axs[1].set_xlim(-80,1350)
	plt.subplots_adjust(left = 0.15,bottom = 0.17,right=0.9)

def plot_MNjeans(files1and2):
        Nmax=0
        fig,axs=plt.subplots(nrows=2, ncols=2,sharex='col')
        plt.subplots_adjust(wspace=0.3, hspace=0)
        colors='fuchsia','k','c'
        labels='8 cells','16 cells','32 cells'
        labels2=r'$\rho_{sink}$=10$^{-10}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-7}$gcm$^{-3}$'
        for j in range(2):
                for i in range((len(files1and2[0]))):
                        N,M,t=txtread(files1and2[j][i])
                        if i==0:
                                  t0=t[0]
                        axs[0,j].plot((t-t0)*code_units.t_cu/(60*60*24*365),N,color=colors[i],label=labels[i])
                        axs[1,j].plot((t-t0)*code_units.t_cu/(60*60*24*365),M*code_units.M_cu/ap.M_sun.cgs.value,color=colors[i],label=labels[i])
                        if N.max()>Nmax:
                                  Nmax=N.max()

                axs[1,j].tick_params(axis="x", labelsize=10,direction="in")
                axs[0,j].tick_params(axis="x", labelsize=10,direction="in")
                axs[0,j].tick_params(axis="y", labelsize=10,direction="in")
                axs[1,j].tick_params(axis="y", labelsize=10,direction="in")
                axs[1,j].text(0.85,0.1,labels2[j],ha='center', va='center', transform=axs[1,j].transAxes,fontsize=10)

        axs[1,0].set_xlabel('t \ [yrs]',fontsize=11)
        axs[1,1].set_xlabel('t \ [yrs]',fontsize=11)
        axs[0,0].set_ylabel(r'N$_{\rm sinks}$',fontsize=11)
        axs[0,1].set_ylabel(r'N$_{\rm sinks}$',fontsize=11)
        axs[1,0].set_ylabel(r'$\sum$ M$_{\rm sink}$ \ [M$_{\odot}]$',fontsize=11)
        axs[1,1].set_ylabel(r'$\sum$ M$_{\rm sink}$ \ [M$_{\odot}]$',fontsize=11)
        #axs[0,0].set_yticks(np.arange(1,Nmax+2,2))
        axs[1,1].legend(fontsize=10,loc='upper left',frameon=False)
        axs[1,0].legend(fontsize=10,loc='upper left',frameon=False)
        #axs[1].set_xlim(-100,1350)
        plt.subplots_adjust(left = 0.1,bottom = 0.17,right=0.95)

def plot_MN_join(files):
        fig = plt.figure()
        ax1=plt.subplot2grid((2, 3), (0, 0), colspan=2)
        ax2=plt.subplot2grid((2, 3), (1, 0), colspan=2,sharex=ax1)
        ax3=plt.subplot2grid((2, 3), (0, 2), rowspan=2)
        
        Nmax=0
        plt.subplots_adjust(wspace=0.4, hspace=0)
        colors='b','g','r','cyan','Purple'
        labels=r'$\rho_{sink}$=10$^{-10}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-9}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-8}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-7}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-6}$gcm$^{-3}$'
        for i in range(len(files)):
                N,M,t=txtread(files[i])
                if i==0:
                        t0=t[0]
                ax1.plot((t-t0)*code_units.t_cu/(60*60*24*365),N,color=colors[i],label=labels[i])
                ax2.plot((t-t0)*code_units.t_cu/(60*60*24*365),M*code_units.M_cu/ap.M_sun.cgs.value,color=colors[i],label=labels[i])
                ax3.plot((t-t0)*code_units.t_cu/(60*60*24*365),N,color=colors[i],label=labels[i])
                if N.max()>Nmax:
                        Nmax=N.max()

        ax2.set_xlabel('t \ [yrs]',fontsize=11)
        ax1.set_ylabel(r'N$_{\rm sinks}$',fontsize=11)
        ax3.set_xlabel('t \ [yrs]',fontsize=11)
        ax3.set_ylabel(r'N$_{\rm sinks}$',fontsize=11)
        ax2.set_ylabel(r'$\sum$ M$_{\rm sink}$ \ [M$_{\odot}]$',fontsize=11)
        ax1.set_yticks(np.arange(1,Nmax+2,2))
        ax2.tick_params(axis="x", labelsize=10,direction="in")
        ax1.tick_params(axis="x", labelsize=10,direction="in")
        ax1.tick_params(axis="y", labelsize=10,direction="in")
        ax2.tick_params(axis="y", labelsize=10,direction="in")
        ax3.tick_params(axis="y", labelsize=10,direction="in")
        ax3.tick_params(axis="x", labelsize=10,direction="in")
        ax3.legend(fontsize=10,loc='upper left',frameon=False,bbox_to_anchor=(0.99, 1.05))
        ax2.set_xlim(-80,1350)
        ax2.set_ylim(-3,70)
        ax3.set_xlim(-50,450)
        ax3.set_ylim(0.5,15)
        plt.subplots_adjust(left = 0.1,bottom = 0.17,right=0.8)
        


def velocity_graph(files):
        colors='Purple','cyan','r','g','b'
        labels=r'$\rho_{sink}$=10$^{-6}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-7}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-8}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-9}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-10}$gcm$^{-3}$'
        for i in range(len(colors)):
                a=aread(files[i])
                x=sum(a.sinkmass*a.sinkx)/sum(a.sinkmass)
                y=sum(a.sinkmass*a.sinky)/sum(a.sinkmass)
                z=sum(a.sinkmass*a.sinkz)/sum(a.sinkmass)
                r=np.sqrt((a.x-x)**2+(a.y-y)**2+(a.z-z)**2)
                mask=np.where(r<1)
                v=np.sqrt(a.vx**2+a.vy**2+a.vz**2) * code_units.v_cu/1e5
                hist, bin_edges = np.histogram(v[mask],weights=a.mass[mask]*code_units.M_cu,bins=np.linspace(0,500,20))
                plt.bar(bin_edges[:-1],np.log10(hist),color=colors[i],width=(bin_edges[1]-bin_edges[0]),label=labels[i])
        plt.xlabel(r'$v \ [kms^{-1}]$',fontsize=20)
        plt.ylabel(r'log$_{10}$(Frequency)',fontsize=20)
        plt.tick_params(axis="x", labelsize=15,direction="in")
        plt.tick_params(axis="y", labelsize=15,direction="in")
        plt.subplots_adjust(left = 0.2,bottom = 0.17,right=0.9)
        plt.tick_params(axis="x", labelsize=15,direction="in")
        plt.tick_params(axis="y", labelsize=15,direction="in")
	#plt.legend(fontsize=12,frameon=False,markerscale=10)
	#I want the legend to read backwards
	
        line1=Line2D([0], [0], color='b', lw=2)
        line2=Line2D([0], [0], color='g', lw=2)
        line3=Line2D([0], [0], color='r', lw=2)
        line4=Line2D([0], [0], color='cyan', lw=2)
        line5=Line2D([0], [0], color='purple', lw=2)
        plt.legend([line1,line2,line3,line4,line5],(r'$\rho_{sink}$=10$^{-10}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-9}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-8}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-7}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-6}$gcm$^{-3}$'),fontsize=12,frameon=False,markerscale=10)



'''||||||||||| radial graphs |||||||||||'''

def vel(files):
	colors='b','g','r','cyan','purple'
	labels=r'$\rho_{sink}$=10$^{-10}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-9}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-8}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-7}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-6}$gcm$^{-3}$'
	fig,ax=plt.subplots(1)
	for i in range(len(files)):
		a=arepo_utils.aread(files[i])
		v=np.sqrt(a.vx**2+a.vy**2+a.vz**2) *code_units.v_cu/1e5
		mid=np.where(a.rho==a.rho.max())
		r=np.sqrt((a.x-a.x[mid])**2+(a.y-a.y[mid])**2+(a.z-a.z[mid])**2) *code_units.d_cu /ap.au.cgs.value
		vs,rs,z=binned_statistic(r,v,bins=10**np.linspace(np.log10(np.sort(r)[1]),np.log10(r.max()),50))
		ax.loglog(rs[:-1],vs,colors[i],label=labels[i])
		ax.legend(fontsize=12,frameon=False)	
		plt.tick_params(axis="x", labelsize=10,direction="in")
		plt.tick_params(axis="y", labelsize=10,direction="in")
		plt.xlabel('R [AU]',fontsize=11)
		plt.ylabel(r'$v$ [kms$^{-1}$]',fontsize=11)

def weighted_average(x,y,weight,bins):
	vals=np.zeros(len(bins)-1)
	for i in range(len(bins)-1):
		mask=np.where((x>=bins[i]) & (x<bins[i+1]))
		if len(mask[0])>0:
			weighted=sum(y[mask]*weight[mask])/sum(weight[mask])
			vals[i]=weighted
		else:
			vals[i]=np.nan
	return vals

def radial(files):
	colors='b','g','r','cyan','purple'
	labels=r'$\rho_{sink}$=10$^{-10}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-9}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-8}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-7}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-6}$gcm$^{-3}$'
	fig,ax=plt.subplots(4,sharex=True)
	plt.subplots_adjust(hspace=0,top=0.95,bottom=0.12,right=0.7,left=0.12)
	for i in range(len(files)):
		a=arepo_utils.aread(files[-1-i])
		v=np.sqrt(a.vx**2+a.vy**2+a.vz**2) *code_units.v_cu/1e5
		midx,midy,midz=a.sinkx[np.where(a.sinkmass==a.sinkmass.max())],a.sinky[np.where(a.sinkmass==a.sinkmass.max())],a.sinkz[np.where(a.sinkmass==a.sinkmass.max())]
		r=np.sqrt((a.x-midx)**2+(a.y-midy)**2+(a.z-midz)**2) *code_units.d_cu /ap.au.cgs.value
	
		rs=10**np.linspace(np.log10(np.sort(r)[1]),np.log10(r.max()),50)
		rhos=weighted_average(r, a.rho*code_units.rho_cu, a.mass/a.rho, rs)
		vs=weighted_average(r, v, a.mass/a.rho, rs)
		Ts=weighted_average(r, a.temp, a.mass/a.rho, rs)
		size=weighted_average(r, ((a.mass/a.rho)**(1/3))*code_units.d_cu, a.mass/a.rho, rs)


		ax[0].loglog(rs[:-1],rhos,colors[-1-i],label=labels[-1-i])
		ax[1].loglog(rs[:-1],vs,colors[-1-i])
		ax[2].semilogx(rs[:-1], Ts,colors[-1-i])
		ax[3].loglog(rs[:-1], size,colors[-1-i])
	for i in range(4):
		ax[i].tick_params(axis="x", labelsize=10,direction="in")
		ax[i].tick_params(axis="y", labelsize=10,direction="in")
	ax[0].set_ylabel(r'$\rho$ [gcm$^{-3}$]',fontsize=10)
	ax[1].set_ylabel(r'v [kms$^{-1}$]',fontsize=10)
	ax[2].set_ylabel('T [K]',fontsize=10)
	ax[3].set_ylabel('L [cm]',fontsize=10)
	ax[3].set_xlabel('R [cm]',fontsize=10)

	line1=Line2D([0], [0], color='b', lw=2)
	line2=Line2D([0], [0], color='g', lw=2)
	line3=Line2D([0], [0], color='r', lw=2)
	line4=Line2D([0], [0], color='cyan', lw=2)
	line5=Line2D([0], [0], color='purple', lw=2)
	ax[0].legend([line1,line2,line3,line4,line5],(r'$\rho_{sink}$=10$^{-10}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-9}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-8}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-7}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-6}$gcm$^{-3}$'),fontsize=10,frameon=False,markerscale=10,loc=(1.01,-0.4))









'''||||||||||| Stromgren sphere function - 1 time use kinda thing ||||||||||||'''



def read_Hirano2014(txtfile):
	r=[]
	M=[]        
	with open(txtfile) as f:
		for line in f.readlines():
			M.append(line.split()[0])
			r.append(line.split()[1])
	return np.asarray(M).astype(float) ,np.asarray(r).astype(float)

def plot_Hirano2014():
	fig,ax=plt.subplots(1)
	M1,r1=read_Hirano2014('P1_curve.txt')
	ax.loglog(M1,r1,'r',label='P1')
	M2,r2=read_Hirano2014('P2_curve.txt')
	ax.loglog(M2,r2,'b',label='P2')
	M3,r3=read_Hirano2014('P3_curve.txt')
	ax.loglog(M3,r3,'k',label='P3')
	ax.tick_params(axis="x", labelsize=10,direction="in",which='both')
	ax.tick_params(axis="y", labelsize=10,direction="in",which='both')
	plt.ylabel(r'R [R$_\odot$]',fontsize=10)	
	plt.xlabel(r'M [M$_\odot$]',fontsize=10)
	plt.legend(fontsize=10,frameon=False)

def stromgren():
	'''Im doing this one in regular units because its monday'''
	fig,ax=plt.subplots(3,sharex=True)
	plt.subplots_adjust(hspace=0)
	colors='b','g','r','cyan','purple'


	#import the data from Hirano2014 
	M1,r1=read_Hirano2014('P1_curve.txt')
	f1=interp1d(M1,r1)

	M2,r2=read_Hirano2014('P2_curve.txt')
	f2=interp1d(M2,r2)

	M3,r3=read_Hirano2014('P3_curve.txt')
	f3=interp1d(M3,r3)

	times=np.array([200,400,600,800])
	alpha=2.6e-13


	Rsphere600=np.array([])
	Rs600=np.array([f3(25),f2(18),f2(13),f1(8),f1(8)]) * ap.R_sun.value
	Ms600=np.array([25,18,13,8,8]) * ap.M_sun.value
	acc600=np.array([6e-3,2e-2,6e-3,1e-2,1e-2]) * ap.M_sun.value /(60*60*24*365)
	L600=ap.G.value * acc600 *Ms600 / Rs600
	T600=(L600/(4*np.pi * Rs600**2*ap.sigma_sb.value))**(1/4)


	Rsphere400=np.array([])
	Rs400=np.array([f3(21),f2(14),f2(9),f1(7),f1(7)]) * ap.R_sun.value
	Ms400=np.array([21,14,9,7,7]) * ap.M_sun.value
	acc400= np.array([2e-2,2e-2,1e-2,4e-3,1e-3])* ap.M_sun.value /(60*60*24*365)
	L400=ap.G.value * acc400 *Ms400 / Rs400
	T400=(L400/(4*np.pi * Rs400**2*ap.sigma_sb.value))**(1/4)
	

	Rsphere200=np.array([])
	Rs200=np.array([f3(16),f2(9),f2(6),f1(5.5),f1(5.5)]) * ap.R_sun.value
	Ms200=np.array([16,9,6,5.5,5.5]) * ap.M_sun.value
	acc200= np.array([3e-2,2e-2,2e-2,5e-3,5e-3])* ap.M_sun.value /(60*60*24*365)
	L200=ap.G.value * acc200 *Ms200 / Rs200
	T200=(L200/(4*np.pi * Rs200**2*ap.sigma_sb.value))**(1/4)

	Rsphere800=np.array([])
	Rs800=np.array([f3(28),f2(19),f2(14),f1(8.5),f1(8.5)]) * ap.R_sun.value
	Ms800=np.array([28,19,14,8.5,8.5]) * ap.M_sun.value
	acc800= np.array([7e-3,4e-3,7e-3,7e-3,7e-3])* ap.M_sun.value /(60*60*24*365)
	L800=ap.G.value * acc800 *Ms800 / Rs800
	T800=(L800/(4*np.pi * Rs800**2*ap.sigma_sb.value))**(1/4)


	n=np.array([1e8,1e9,1e10,1e11,1e12])*code_units.rho_cu/ap.m_p.cgs.value   #allowed to be in cgs because alpha is 


	for i in range(5):
		v=10**np.linspace(np.log10(3.28e15),20,100) #visible to 
		plank600=2*ap.h.value*v**3/ap.c.value**2 *1/(np.exp(ap.h.value*v/(ap.k_B.value*T600[i]))-1)
		dv=np.zeros_like(v)
		dv[1:]=v[1:]-v[:-1]
		N600=sum(4*np.pi * plank600*dv /(ap.h.value*v) * 4*np.pi*Rs600[i]**2)
		Rsphere600=np.append(Rsphere600,(N600 * 3/(4*np.pi*n[i]**2*alpha))**(1/3) /100) #no cgs allowed

		plank400=2*ap.h.value*v**3/ap.c.value**2 *1/(np.exp(ap.h.value*v/(ap.k_B.value*T400[i]))-1)		
		N400=sum(4*np.pi * plank400*dv /(ap.h.value*v) * 4*np.pi*Rs400[i]**2)
		Rsphere400=np.append(Rsphere400,(N400 * 3/(4*np.pi*n[i]**2*alpha))**(1/3) /100)

		plank200=2*ap.h.value*v**3/ap.c.value**2 *1/(np.exp(ap.h.value*v/(ap.k_B.value*T200[i]))-1)
		N200=sum(4*np.pi * plank200*dv /(ap.h.value*v) * 4*np.pi*Rs200[i]**2)
		Rsphere200=np.append(Rsphere200,(N200 * 3/(4*np.pi*n[i]**2*alpha))**(1/3) /100)

		plank800=2*ap.h.value*v**3/ap.c.value**2 *1/(np.exp(ap.h.value*v/(ap.k_B.value*T800[i]))-1)
		N800=sum(4*np.pi * plank800*dv /(ap.h.value*v) * 4*np.pi*Rs800[i]**2)
		Rsphere800=np.append(Rsphere800,(N800 * 3/(4*np.pi*n[i]**2*alpha))**(1/3) /100)

	for i in range(5):
		linestyle='solid'
		if i==4:
			linestyle='--'
		ax[0].plot(times,np.array([Rs200[i],Rs400[i],Rs600[i],Rs800[i]])/ap.R_sun.value,color=colors[i],linestyle=linestyle,marker='o')

		#ax[1].plot(times,np.array([L200[i],L400[i],L600[i],L800[i]])/ap.L_sun.value,color=colors[i])

		ax[1].plot(times,np.array([T200[i],T400[i],T600[i],T800[i]]),color=colors[i],linestyle=linestyle,marker='o')

		ax[2].semilogy(times,np.array([Rsphere200[i],Rsphere400[i],Rsphere600[i],Rsphere800[i]])/ap.R_sun.cgs.value,color=colors[i],linestyle=linestyle,marker='o')
		if i==4:
			ax[0].plot(times[:2],np.array([Rs200[i],Rs400[i]])/ap.R_sun.value,color=colors[i])
			ax[1].plot(times[:2],np.array([T200[i],T400[i]]),color=colors[i])
			ax[2].semilogy(times[:2],np.array([Rsphere200[i],Rsphere400[i]])/ap.R_sun.cgs.value,color=colors[i])
	for i in range(3):
		ax[i].tick_params(axis="x", labelsize=10,direction="in",which='both')
		ax[i].tick_params(axis="y", labelsize=10,direction="in",which='both')
	ax[2].set_xlabel('t [yrs]',fontsize=10)
	ax[0].set_ylabel(r'R$_\star$ [R$_{\odot}$]',fontsize=10)
	ax[1].set_ylabel('T [K]',fontsize=10)
	ax[2].set_ylabel(r'R$_S$ [R$_{\odot}$]',fontsize=10)
			



'''||||||||||| functions for looking at the largest sink vs time ||||||||||||'''

def write_largest_sink(dirname,start,end,interval,name):
	if os.path.isfile(name):
		print("file exists, appending")
		f=open(name, "a+")
	else:
		print("creating new file")
		f = open(name, "x")
	num=int((end-start)/interval)
	counter=0
	for i in range(num):
		text_trap = io.StringIO() #prevent massive text output from snapshot reads
		sys.stdout = text_trap
		n=snapname(start,i,interval)
		a=arepo_utils.aread(dirname+'snapshot_'+n)
		sys.stdout = sys.__stdout__
	
			
		if a.npart[-1]>0:
				counter+=1
				M=a.sinkmass.max() * code_units.M_cu / ap.M_sun.cgs.value
				t=a.time * code_units.t_cu / (60*60*24*365)
				if counter==1:
					f.write(str(M)+' '+str(0)+' '+str(t) + '\n')
					massold=M
					timeold=t
				else:
					accrate= (M-massold) / (t-timeold)
					f.write(str(M)+' '+str(accrate)+' '+str(t) + '\n')
					massold=M
					timeold=t
		print(str(n)+' :done')
	f.close()
	
	
	
def read_largest_sink(txtfile):
	M=[]
	acc=[]
	t=[]
	with open(txtfile) as f:
		for line in f.readlines():
			M.append(line.split()[0])
			acc.append(line.split()[1])
			t.append(line.split()[2])
	return np.asarray(M).astype(float) ,np.asarray(acc).astype(float), np.asarray(t).astype(float)



def plot_largest_sink(files):

	colors='b','g','r','cyan','Purple'
	labels=r'$\rho_{sink}$=10$^{-10}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-9}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-8}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-7}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-6}$gcm$^{-3}$'

	fig,ax =plt.subplots(2,sharex=True)
	plt.subplots_adjust(hspace=0)
	for i in range(len(files)):
		M,acc,t=read_largest_sink(files[i])
		#if i==0:	
		t0=t[0]
		if i==0:
			bins=np.linspace(t0,t.max(),200)
		acc,T,z=binned_statistic(t,acc,bins=bins)
		interval=5
		#if i==4:
		#	interval=10
		#mask=np.where(acc[1::interval]>0)
		ax[0].semilogy(T[:-1]-t0,acc,colors[i])
		ax[1].semilogy(t-t0,M,colors[i],label=labels[i])

	plt.subplots_adjust(left = 0.2,bottom = 0.17,right=0.9)
	ax[0].tick_params(axis="x", labelsize=10,direction="in")
	ax[0].tick_params(axis="y", labelsize=10,direction="in")
	ax[1].tick_params(axis="x", labelsize=10,direction="in")
	ax[1].tick_params(axis="x", labelsize=10,direction="in")
	ax[1].set_xlabel('t [yrs]',fontsize=10)
	ax[0].set_ylabel(r'$\dot {\rm M}_{\rm largest}$ [M$_\odot$ yr$^{-1}$]',fontsize=10)
	ax[1].set_ylabel(r'M$_{\rm largest}$ [M$_\odot$]',fontsize=10)
	ax[1].legend(frameon=False,loc='lower right',fontsize=10)
	#plt.axhline(y=8,color='k')
	return fig,ax


def MMdot(files):
	for i in range(5):
		colors='b','g','r','cyan','Purple'
		labels=r'$\rho_{sink}$=10$^{-10}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-9}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-8}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-7}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-6}$gcm$^{-3}$'
		M,acc,t=read_largest_sink(files[i])
		acc,M,b=binned_statistic(M,acc,bins=10**np.linspace(0,2,20))
		mask=np.where(acc>0)
		plt.loglog(M[:-1][mask],acc[mask],c=colors[i])
	plt.axhline(y=0.004,color='k',linestyle='--')
	plt.axhline(y=0.04,color='k',linestyle='--')
	plt.tick_params(axis="x", labelsize=10,direction="in")
	plt.tick_params(axis="y", labelsize=10,direction="in")
	plt.ylabel(r'$\rm \dot M_{\rm largest}}$ [M$_\odot$yr$^{-1}$]',fontsize=10)
	plt.xlabel(r'M [M$_\odot$]',fontsize=10)
