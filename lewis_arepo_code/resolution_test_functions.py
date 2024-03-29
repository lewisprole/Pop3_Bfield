
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
import read_sink_info
from scipy.integrate import quad


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
        
'''ejection checker'''

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
	print('no eject: '+str(ejects))
	print(masses)
	return args.astype(int)

def ejection_fraction(files):
	bins=np.array([0,0.075,0.8,1000])
	mask1=np.array([])
	mask_ej1=np.array([])
	mask2=np.array([])
	mask_ej2=np.array([])
	mask3=np.array([])
	mask_ej3=np.array([])
	fracs=np.array([])
	for i in range(len(files)):
		a=arepo_utils.aread(files[i])
		fracs=np.array([])
		args=energy_check(a,0.03)
			
		mask1=np.append(mask1, (np.where((a.sinkmass>bins[0]) & (a.sinkmass<=bins[1])))[0] )
		mask_ej1=np.append(mask_ej1, (np.where((a.sinkmass[args]>bins[0]) & (a.sinkmass[args]<=bins[1])))[0])

		mask2=np.append(mask2, (np.where((a.sinkmass>bins[1]) & (a.sinkmass<=bins[2])))[0] )
		mask_ej2=np.append(mask_ej2, (np.where((a.sinkmass[args]>bins[1]) & (a.sinkmass[args]<=bins[2])))[0])

		mask3=np.append(mask3, (np.where((a.sinkmass>bins[2]) & (a.sinkmass<=bins[3])))[0] )
		mask_ej3=np.append(mask_ej3, (np.where((a.sinkmass[args]>bins[2]) & (a.sinkmass[args]<=bins[3])))[0])

	fracs=np.append(fracs, len(mask_ej1)/len(mask1))
	fracs=np.append(fracs, len(mask_ej2)/len(mask2))
	fracs=np.append(fracs, len(mask_ej3)/len(mask3))
	return fracs
	
def ejction_time(dirname,start,end,interval,min_eject_length):
	Neject=np.array([])
	t=np.array([])
	name=dirname+'/ejections.txt'
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
		args=energy_check(a,min_eject_length)
		Neject=np.append(Neject,len(args))
		t=np.append(t,a.time)
		f.write(str(a.time) + ' ' + str(len(args)) + ' \n')
		print(str(n)+' :done')
	f.close()

def read_ejection(filename):
	t=np.array([])
	Neject=np.array([])
	with open(filename) as f:
		for line in f.readlines():
			t=np.append(t,float(line.split()[0]))
			Neject=np.append(Neject,float(line.split()[1]))
	return t,Neject


def ejection_graph():
	fig,ax=plt.subplots(1)
	dirs='/scratch/c.c1521474/resolution_test/merge/1e12/','/scratch/c.c1521474/resolution_test/seed4/1e12/','/scratch/c.c1521474/resolution_test/seed5/1e12/'
	for i in range(3):
		t,N,M=read_sink_info.Nsinks(dirs[i]+'/sink_particle_info/')
		T,Neject=read_ejection(dirs[i]+'/ejections.txt')
		Tnew=np.array([T[0]])
		Nnew=np.array([Neject[0]])
		for j in range(len(T)-1):
			if Neject[j+1]>Neject[j]:
				Tnew=np.append(Tnew,T[j+1])
				Tnew=np.append(Tnew,T[j+1])
				Nnew=np.append(Nnew,Neject[j])
				Nnew=np.append(Nnew,Neject[j+1])
			if j==len(T)-2:
				Tnew=np.append(Tnew,T[j+1])
				Nnew=np.append(Nnew,Neject[j+1])
		ax.plot((Tnew-t[0])*code_units.t_cu/(60*60*24*365),Nnew,label=('A','B','C')[i])
	plt.subplots_adjust(left=0.2,right=0.8,top=0.8,bottom=0.2)
	plt.subplots_adjust(hspace=0)	
	ax.tick_params(axis="y", labelsize=10,direction="in",which='both')
	ax.tick_params(axis="x", labelsize=10,direction="in",which='both')
	ax.set_xlabel('t [yr]')
	ax.set_ylabel(r'N$_{\rm eject}$')
	ax.legend(frameon=False,fontsize=8)

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

	times=np.array([200,400,600,800,1000])
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

	Rsphere1000=np.array([])
	Rs1000=np.array([f3(28.6),f2(19.5),f2(15.2)]) * ap.R_sun.value
	Ms1000=np.array([28.6,19.5,15.2]) * ap.M_sun.value
	acc1000= 5e-3* ap.M_sun.value /(60*60*24*365)
	L1000=ap.G.value * acc1000 *Ms1000 / Rs1000
	T1000=(L1000/(4*np.pi * Rs1000**2*ap.sigma_sb.value))**(1/4)


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
		if i<3:
			plank1000=2*ap.h.value*v**3/ap.c.value**2 *1/(np.exp(ap.h.value*v/(ap.k_B.value*T1000[i]))-1)
			N1000=sum(4*np.pi * plank1000*dv /(ap.h.value*v) * 4*np.pi*Rs1000[i]**2)
			Rsphere1000=np.append(Rsphere1000,(N1000 * 3/(4*np.pi*n[i]**2*alpha))**(1/3) /100)
	for i in range(2):
		T1000=np.append(T1000,np.nan)
		Rs1000=np.append(Rs1000,np.nan)
		Rsphere1000=np.append(Rsphere1000,np.nan)

	for i in range(5):
		linestyle='solid'
		if i==4:
			linestyle='--'
		ax[0].plot(times,np.array([Rs200[i],Rs400[i],Rs600[i],Rs800[i],Rs1000[i]])/ap.R_sun.value,color=colors[i],linestyle=linestyle,marker='o')

		#ax[1].plot(times,np.array([L200[i],L400[i],L600[i],L800[i]])/ap.L_sun.value,color=colors[i])

		ax[1].plot(times,np.array([T200[i],T400[i],T600[i],T800[i],T1000[i]]),color=colors[i],linestyle=linestyle,marker='o')

		ax[2].semilogy(times,np.array([Rsphere200[i],Rsphere400[i],Rsphere600[i],Rsphere800[i],Rsphere1000[i]])/ap.R_sun.cgs.value,color=colors[i],linestyle=linestyle,marker='o')
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

	ax[2].set_ylim(5e-14,2e-5)
	ax[1].set_ylim(1500,8500)
	ax[0].set_ylim(0,700)
			



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
		#if i==0:
		#	bins=np.linspace(t0,t.max(),200)
		#acc,T,z=binned_statistic(t,acc,bins=bins)
		#interval=5
		T=t
		if i==4:
			T=T[0::3][np.where(acc[0::3]>0)]
			acc=acc[0::3][np.where(acc[0::3]>0)]
		#mask=np.where(acc[1::interval]>0)
		ax[0].semilogy(T-t0,acc,colors[i])
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
		acc,M,b=binned_statistic(M,acc,bins=10**np.linspace(0,2,100))
		mask=np.where(acc>0)
		plt.loglog(M[:-1][mask],acc[mask],c=colors[i])
	plt.axhline(y=0.004,color='k',linestyle='--')
	plt.axhline(y=0.04,color='k',linestyle='--')
	plt.tick_params(axis="x", labelsize=10,direction="in")
	plt.tick_params(axis="y", labelsize=10,direction="in")
	plt.ylabel(r'$\rm \dot M_{\rm largest}}$ [M$_\odot$yr$^{-1}$]',fontsize=10)
	plt.xlabel(r'M [M$_\odot$]',fontsize=10)





def stromgren_spheres(files):

	fig,ax=plt.subplots(5,sharex=True)
	plt.subplots_adjust(hspace=0,left=0.15,right=0.7,top=0.95,bottom=0.08)	
	colors='b','g','r','cyan','purple'
	labels=r'$\rho_{sink}$=10$^{-10}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-9}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-8}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-7}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-6}$gcm$^{-3}$'
	
	M1,r1=read_Hirano2014('P1_curve.txt')
	f1=interp1d(M1,r1,fill_value='extrapolate')
	M2,r2=read_Hirano2014('P2_curve.txt')
	f2=interp1d(M2,r2,fill_value='extrapolate')
	M3,r3=read_Hirano2014('P3_curve.txt')
	f3=interp1d(M3,r3,fill_value='extrapolate')

	times=np.array([200,400,600,800,1000])
	alpha=2.6e-13
	n=np.array([1e8,1e9,1e10,1e11,1e12])*code_units.rho_cu/ap.m_p.cgs.value #ok this one needs to be cgs becasue of alpha 
	for i in range(len(files)):
		M,acc,t=read_largest_sink(files[i])
		if i==len(files)-1:
			M,acc,t=M[0::3],acc[0::3],t[0::3]
		if i==0:
			rs=f3(M) *ap.R_sun.value
		if (i>0) & (i<3):
			rs=f2(M) *ap.R_sun.value
		if i>2:
			rs=f1(M) *ap.R_sun.value
		M=M*ap.M_sun.value
		acc=acc*ap.M_sun.value/(60*60*24*365)
		L=ap.G.value * acc *M / rs
		T=(L/(4*np.pi * rs**2*ap.sigma_sb.value))**(1/4)	


		v=10**np.linspace(np.log10(3.28e15),20,100) #Lynman to xray
		dv=np.zeros_like(v)
		dv[1:]=v[1:]-v[:-1]
		Rsphere=np.array([])
		for j in range(len(T)):
			plank=2*ap.h.value*v**3/ap.c.value**2 *1/(np.exp(ap.h.value*v/(ap.k_B.value*T[j]))-1)
			N=sum(4*np.pi * plank*dv /(ap.h.value*v) * 4*np.pi*rs[j]**2)
			Rsphere=np.append(Rsphere,(N * 3/(4*np.pi*n[i]**2*alpha))**(1/3) /100) #no cgs allowed
		mask=np.where(acc>0)
		ax[0].semilogy(t[mask]-t[0],acc[mask]/(ap.M_sun.value/(60*60*24*365)),color=colors[i],label=labels[i])
		ax[1].semilogy(t-t[0],M/ap.M_sun.value,color=colors[i])
		ax[2].plot(t-t[0],rs/ap.R_sun.value,color=colors[i])
		ax[3].plot(t-t[0],T,color=colors[i])
		ax[4].semilogy(t-t[0],Rsphere/ap.R_sun.cgs.value,colors[i])
	ax[0].legend(fontsize=10,frameon=False,markerscale=10,loc=(1.01,-0.2))
	ax[0].set_ylabel(r'$\rm \dot M$ [M$_\odot$yr$^{-1}$]',fontsize=10)
	ax[1].set_ylabel(r'M [M$_\odot$]',fontsize=10)
	ax[2].set_ylabel(r'R$_\star$ [R$_\odot$]',fontsize=10)
	ax[3].set_ylabel('T [K]',fontsize=10)
	ax[4].set_ylabel(r'R$_S$ [R$_\odot$]',fontsize=10)
	ax[4].set_xlabel('t [yrs]')


''' MHD functions '''

def hydro_checkcrop(files_crop,files_big):
	fig,ax=plt.subplots(1,3,sharey=True)
	plt.subplots_adjust(wspace=0)
	colors='b','g','r'
	titles=r'10$^{-10}$gcm$^{-3}$',r'10$^{-7}$gcm$^{-3}$',r'10$^{-8}$gcm$^{-3}$'
	ax[0].set_ylabel(r'N$_{\rm sinks}$',fontsize=11)
	for i in range(len(files_crop)):
		N,M,t=txtread(files_crop[i])
		ax[i].plot((t-t[0])*code_units.t_cu/(60*60*24*365),N,color=colors[i],label='crop')
		N,M,t=txtread(files_big[i])
		ax[i].plot((t-t[0])*code_units.t_cu/(60*60*24*365),N,linestyle='--',color=colors[i],label='full')
		ax[i].set_xlabel('t [yrs]',fontsize=11)
		ax[i].tick_params(axis="x", labelsize=10,direction="in",which='both')
		ax[i].tick_params(axis="y", labelsize=10,direction="in",which='both')
		ax[i].legend(fontsize=10,frameon=False,loc='upper left')
		ax[i].set_title(titles[i],fontsize=10)
		


def MHD_compare(files_noMHD,files_uniform,files_MHD):
	fig,ax=plt.subplots(1,3,sharey=True)
	plt.subplots_adjust(wspace=0)
	colors='g','r','cyan'
	ax[0].set_ylabel(r'N$_{\rm sinks}$',fontsize=11)
	for i in range(len(files_noMHD)):
		N,M,t=txtread(files_noMHD[i])
		ax[i].plot((t-t[0])*code_units.t_cu/(60*60*24*365),N,color='aquamarine',label='hydro')
		N,M,t=txtread(files_MHD[i])
		ax[i].plot((t-t[0])*code_units.t_cu/(60*60*24*365),N,color='royalblue',label=r'$\propto$k$^{3/2}$')
		N,M,t=txtread(files_uniform[i])
		ax[i].plot((t-t[0])*code_units.t_cu/(60*60*24*365),N,color='Crimson',label='uniform field')
		
		ax[i].set_xlabel('t [yrs]',fontsize=11)
		ax[i].tick_params(axis="x", labelsize=10,direction="in",which='both')
		ax[i].tick_params(axis="y", labelsize=10,direction="in",which='both')
		if i==0:
			ax[i].legend(fontsize=10,frameon=False,loc='upper left')

def Bfield_plot(files,sink):
	fig,ax=plt.subplots(2,sharex=True)
	plt.subplots_adjust(hspace=0)
	colors='cyan','r','g'
	labels=r'$\rho_{sink}$=10$^{-9}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-8}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-7}$gcm$^{-3}$'
	for i in range(len(files)):	
		a=arepo_utils.aread(files[i])
		B=np.sqrt(a.bfield[:,0]**2+ a.bfield[:,1]**2 + a.bfield[:,2]**2)
		#if a.npart[-1]>0:
			#midx=a.sinkx[np.where(a.sinkmass==a.sinkmass.max())]
			#midy=a.sinky[np.where(a.sinkmass==a.sinkmass.max())]
			#midz=a.sinkz[np.where(a.sinkmass==a.sinkmass.max())]
			#r=np.sqrt((a.x-midx)**2+(a.y-midy)**2+(a.z-midz)**2)
		#else:
			#mid = np.where(a.rho==a.rho.max())	
			#r=np.sqrt((a.x-a.x[mid])**2+(a.y-a.y[mid])**2+(a.z-a.z[mid])**2)
		B,rho,y=binned_statistic(a.rho*code_units.rho_cu,B*code_units.B_cu,bins=10**np.linspace(np.log10(a.rho.min()*code_units.rho_cu),np.log10(a.rho.max()*code_units.rho_cu),50),statistic='median')
		ax[0].loglog(rho[:-1],B*code_units.B_cu,color=colors[i])	
		ax[1].loglog(rho[:-1],B*code_units.B_cu/rho[:-1]**(2/3),color=colors[i])
		#B,r,y=binned_statistic(r,B,bins=10**np.linspace(np.log10(np.sort(r)[1]),np.log10(r.max()),50))
		#ax.loglog(r[:-1]*code_units.d_cu,B*code_units.B_cu)
	ax[1].loglog(rho[:-1],rho[:-1]**0.05*2e3,'k')
	ax[0].set_ylabel('B [G]')
	ax[1].set_ylabel(r'B/$\rho^{2/3}$')
	ax[1].set_xlabel(r'$\rho$ [gcm$^{-3}$]')


def angular_velocity(a_sink,a_initial):
	if len(a_sink.sinkx)==1:
		x=a_sink.sinkx-a_initial.x
		y=a_sink.sinky-a_initial.y
		z=a_sink.sinkz-a_initial.z
	else:
		midx=sum(a_sink.sinkx*a_sink.sinkmass)/sum(a_sink.sinkmass)
		midy=sum(a_sink.sinky*a_sink.sinkmass)/sum(a_sink.sinkmass)
		midz=sum(a_sink.sinkz*a_sink.sinkmass)/sum(a_sink.sinkmass)
		x=midx-a_initial.x
		y=midy-a_initial.y
		z=midz-a_initial.z

	r=np.sqrt(x**2+y**2+z**2)
	crossx=a_initial.vy*z - a_initial.vz*y
	crossy=-(a_initial.vx*z - a_initial.vz*x)
	crossz=a_initial.vx*y - a_initial.vy*x
	cross=np.sqrt(crossx**2+crossy**2+crossz**2)
	vrot=cross/r
	v=np.sqrt(a_initial.vx**2+a_initial.vy**2+a_initial.vz**2)
	vs,rs,z=binned_statistic(r,vrot/v,bins=10**np.linspace(np.log10(np.sort(r)[1]),np.log10(r.max()),50))
	return vs,rs

def angular_graph(sink,initial):
	sink='/scratch/c.c1521474/resolution_test/merge/1e12/snapshot_039','/scratch/c.c1521474/resolution_test/seed4/1e12/snapshot_036','/scratch/c.c1521474/resolution_test/seed5/1e12/snapshot_021'
	initial='/scratch/c.c1521474/resolution_test/merge/1e12/snapshot_000','/scratch/c.c1521474/resolution_test/seed4/1e12/snapshot_000','/scratch/c.c1521474/resolution_test/seed5/1e12/snapshot_000'
	fig,ax=plt.subplots(1)
	ax.tick_params(axis="x", labelsize=10,direction="in",which='both')
	ax.tick_params(axis="y", labelsize=10,direction="in",which='both')
	ax.set_ylabel(r'v$_\theta$ / v',fontsize=10)
	ax.set_xlabel('R [pc]',fontsize=10)
	for i in range(3):
		a_sink=arepo_utils.aread(sink[i])
		a_initial=arepo_utils.aread(initial[i])
		vs,rs=angular_velocity(a_sink,a_initial)
		ax.plot(rs[:-1]*code_units.d_cu/ap.pc.cgs.value,vs,label=('A','B','C')[i])
	ax.legend(fontsize=10,loc='lower right',frameon=False)
	plt.subplots_adjust(left=0.2,right=0.8,top=0.8,bottom=0.2)
	ax.set_yscale('log')

def Nphotons(v,T,R):
	plank=2*ap.h.value*v**3/ap.c.value**2 *1/(np.exp(ap.h.value*v/(ap.k_B.value*T))-1)
	return np.pi * plank /(ap.h.value*v) * 4*np.pi*R**2


def radial_stromgren(snapA,dirname):
	'''want to estimate the stromgren sphere radius using different initial spheres to use their average density'''
	a=arepo_utils.aread(snapA)
	#b=arepo_utils.aread(snapB)
	#ACC=np.zeros_like(a.sinkx)
	#M=a.sinkmass
	#for i in range(len(a.sinkx)):
	#	if a.idsink[i] in b.idsink:
	#		arg=np.where(b.idsink==a.idsink[i])
	#		ACC[i]=(a.sinkmass[i]-b.sinkmass[arg])/(a.time-b.time)
	#ACC=ACC*ap.M_sun.value/code_units.t_cu
	#M=M*ap.M_sun.value


	M,ACC,time,x,y,z=read_sink_info.all_mass(dirname)
	current=np.where(abs(time-a.time)==abs(time-a.time).min())[0][0]
	ACC=ACC[current,:]*ap.M_sun.value /(60*60*24*365)
	M=M[current,:]*ap.M_sun.value
	x=x[current,:]
	y=y[current,:]
	z=z[current,:]

	R=26*ap.R_sun.value*(M/ap.M_sun.value)**(0.27) * (ACC *60*60*24*365/ap.M_sun.value/1e-3)**(0.41)

	#accretion luminosity 
	Lacc=M*ACC*ap.G.value/ R

	#Temperature from Stefan-Boltzman
	T=(Lacc/(4*np.pi*R**2) /ap.sigma_sb.value)**(1/4)

	#Plank function
	Nion=np.zeros_like(M)
	Nion_all=0
	for i in range(len(M)):
		if T[i]>0:
			
			v_peak=5.88e10*T[i]
			if v_peak>3.28e15:
				point=v_peak
			if v_peak<=3.28e15:
        			point=3.28e15
			N=quad(Nphotons,3.28e15,1e25,args=(T[i],R[i]),points=np.array([point,point*2,point*5,point*10]))
			Nion[i]=N[0]
			Nion_all+=N[0]
	
	for i in range(len(M)):
		if Nion[i]>0:
			print(i)
			rs=np.sqrt((a.x-x[i])**2+(a.y-y[i])**2+(a.z-z[i])**2)
			r_bins=10**np.linspace(np.log10(np.sort(rs)[1]),np.log10(rs.max()),100)
			R_str=np.zeros_like(r_bins)
			for j in range(len(r_bins)-1):
				mask=np.where(rs<r_bins[j+1])
				rho=np.mean(a.rho[mask]) *code_units.rho_cu * 1000 /ap.m_p.value #base unit
				alpha=2.6e-19
				R_str[j]=(Nion[i] * 3/(4*np.pi*rho**2*alpha))**(1/3)
			print(R_str)
			plt.loglog(r_bins*code_units.d_cu/ap.R_sun.cgs.value,R_str/ap.R_sun.value)
	
	argmax=np.where(M==M.max())[0]
	rs=np.sqrt((a.x-x[argmax])**2+(a.y-y[argmax])**2+(a.z-z[argmax])**2)
	r_bins=10**np.linspace(np.log10(np.sort(rs)[1]),np.log10(rs.max()),100)
	R_str=np.zeros_like(r_bins)
	for j in range(len(r_bins)-1):
		mask=np.where(rs<r_bins[j+1])
		rho=np.median(a.rho[mask]) *code_units.rho_cu * 1000 /ap.m_p.value #base unit
		R_str[j]=(Nion_all * 3/(4*np.pi*rho**2*alpha))**(1/3)
		plt.loglog(r_bins*code_units.d_cu/ap.R_sun.cgs.value,R_str/ap.R_sun.value,'k')

	plt.ylabel(r'R$_{\rm Str}$')
	plt.xlabel(r'R$_{\rm sphere}$')
