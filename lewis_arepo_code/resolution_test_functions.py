
import numpy as np
import matplotlib.pyplot as plt 
import arepo_utils
import io 
import sys
import code_units
import astropy.constants as ap
import os


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
	fig,axs=plt.subplots(2,sharex=True)
	plt.subplots_adjust(wspace=0, hspace=0)
	colors='b','g','r','cyan'
	labels=r'$\rho_{sink}$=10$^{-10}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-9}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-8}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-7}$gcm$^{-3}$'
	for i in range(len(files)):
		N,M,t=txtread(files[i])
		if i==0:
			t0=t[0]
		axs[0].plot((t-t0)*code_units.t_cu/(60*60*24*365),N,color=colors[i])
		axs[1].plot((t-t0)*code_units.t_cu/(60*60*24*365),M*code_units.M_cu/ap.M_sun.cgs.value,color=colors[i],label=labels[i])

	axs[1].set_xlabel(r'$t \ [yrs]$',fontsize=20)
	axs[0].set_ylabel(r'$N_{sinks}$',fontsize=20)
	axs[1].set_ylabel(r'$\sum M_{sink} \ [M_{\odot}]$',fontsize=20)
	axs[1].tick_params(axis="x", labelsize=15)
	axs[0].tick_params(axis="y", labelsize=15)
	axs[1].tick_params(axis="y", labelsize=15)
	axs[1].legend(fontsize=12)

	


		

