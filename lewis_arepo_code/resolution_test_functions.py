
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
	colors='b','g','r','cyan'
	labels=r'$\rho_{sink}$=10$^{-10}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-9}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-8}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-7}$gcm$^{-3}$'
	for i in range(len(files)):
		N,M,t=txtread(files[i])
		if i==0:
			t0=t[0]
		axs[0].plot((t-t0)*code_units.t_cu/(60*60*24*365),N,color=colors[i],label=labels[i])
		axs[1].plot((t-t0)*code_units.t_cu/(60*60*24*365),M*code_units.M_cu/ap.M_sun.cgs.value,color=colors[i],label=labels[i])
		if N.max()>Nmax:
			Nmax=N.max()

	axs[1].set_xlabel(r'$t \ [yrs]$',fontsize=20)
	axs[0].set_ylabel(r'$N_{sinks}$',fontsize=20)
	axs[1].set_ylabel(r'$\sum M_{sink} \ [M_{\odot}]$',fontsize=20)
	axs[0].set_yticks(np.arange(1,Nmax+2,2))
	axs[1].tick_params(axis="x", labelsize=15,direction="in")
	axs[0].tick_params(axis="x", labelsize=15,direction="in")
	axs[0].tick_params(axis="y", labelsize=15,direction="in")
	axs[1].tick_params(axis="y", labelsize=15,direction="in")
	axs[0].legend(fontsize=12,loc='upper left',frameon=False,bbox_to_anchor=(0.99, 1))
	axs[1].set_xlim(-80,1350)
	plt.subplots_adjust(left = 0.15,bottom = 0.17,right=0.7)

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
		t0=t[0]
		acc,T,z=binned_statistic(t,acc,bins=20)
		ax[0].semilogy(T[1:]-t0,acc,colors[i])
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

