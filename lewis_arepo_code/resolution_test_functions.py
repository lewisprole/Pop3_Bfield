
import numpy as np
import matplotlib.pyplot as plt 
import arepo_utils
import io 
import sys
import code_units
import astropy.constants as ap
import os
from matplotlib.lines import Line2D


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

def plot_MNjeans(files):
        Nmax=0
        fig,axs=plt.subplots(2,sharex=True)
        plt.subplots_adjust(wspace=0, hspace=0)
        colors='fuchsia','k','c'
        labels='8 cells','16 cells','32 cells'
        for i in range(len(files)):
                N,M,t=txtread(files[i])
                if i==0:
                        t0=t[0]
                axs[0].plot((t-t0)*code_units.t_cu/(60*60*24*365),N,color=colors[i])
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
        axs[1].legend(fontsize=12,loc='upper left',frameon=False)
        #axs[1].set_xlim(-100,1350)
	

def velocity_graph(files):
        colors='cyan','r','g','b'
        labels=r'$\rho_{sink}$=10$^{-7}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-8}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-9}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-10}$gcm$^{-3}$'
	for i in range(len(colors)):
		a=aread(files[i])
		x=sum(a.sinkmass*a.sinkx)/sum(a.sinkmass)
		y=sum(a.sinkmass*a.sinky)/sum(a.sinkmass)
		z=sum(a.sinkmass*a.sinkz)/sum(a.sinkmass)
		r=np.sqrt((a.x-x)**2+(a.y-y)**2+(a.z-z)**2)
		mask=np.where(r<1)
		v=np.sqrt(a.vx**2+a.vy**2+a.vz**2) * v_cu/1e5
		hist, bin_edges = np.histogram(v[mask],weights=d_cu**3*(a.mass/a.rho)[mask],bins=np.linspace(0,300,50))
		plt.semilogy(bin_edges[:-1],hist,c=colors[i],label=labels[i])
	plt.xlabel(r'$v \ [kms^{-1}]$',fontsize=20)
	plt.ylabel('Volume weighted frequency',fontsize=20)
	plt.tick_params(axis="x", labelsize=15,direction="in")
	plt.tick_params(axis="y", labelsize=15,direction="in")
	plt.subplots_adjust(left = 0.2,bottom = 0.17,right=0.9)
	#plt.legend(fontsize=12,frameon=False,markerscale=10)
	#I want the legend to read backwards
	
	line1=Line2D([0], [0], color='b', lw=2)
	line2=Line2D([0], [0], color='g', lw=2)
	line3=Line2D([0], [0], color='r', lw=2)
	line4=Line2D([0], [0], color='cyan', lw=2)
	line5=Line2D([0], [0], color='purple', lw=2)
	plt.legend([line1,line2,line3,line4],(r'$\rho_{sink}$=10$^{-10}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-9}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-8}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-7}$gcm$^{-3}$'),fontsize=12,frameon=False,markerscale=10)





