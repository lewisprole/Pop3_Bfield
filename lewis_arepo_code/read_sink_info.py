import numpy as np 
import matplotlib.pyplot as plt 
import struct
import code_units
import astropy.constants as ap 
import os
from scipy.interpolate import interp1d
plt.ion()

def allfiles(dirname):
	'''takes all files in directory and orders them based on order of ending 3 number string'''
	files=os.listdir(dirname)
	files_new=np.array([])
	nos=np.array([])
	for i in range(len(files)):
		no=files[i][-3:]
		if no[0]==0: #get rid of zeros e.g. '001' -> '1'
			no=no[1:]
			if no[0]==0:
				no=no[1:]
		no=int(no)
		nos=np.append(nos,no)
	args=np.argsort(nos) #rearrange list 
	for i in range(len(args)):
		files_new=np.append(files_new,files[args[i]])
	return files_new



def Nsinks(dirname):
	'''cycle starts wih double: time, int: Nsinks, then...
	doubles: 3pos, 3vel, 3accel, mass, FormationMass, FormationTime, AccretionRate, TimeOld =14
	long long ID, ints: HomeTask, Index, FormationOrder
	so sink mass starts at (8+4 + 8*9) bits, then every N(8*14 + 8 + 4*3 + 4) where N is the Number of sinks cycled through so far'''
	#first get all the sink_particle_info files into the right order
	files = allfiles(dirname)
	Nsinks_old =0 
	T=np.array([])
	N=np.array([])
	M=np.array([])
	eject=np.array([])
	for i in range(len(files)):
		with open(dirname+files[i], mode='rb') as file:
			data = file.read()
		bits_sink= 8*14 + 8 + 4*3 + 4
		n=0
		cycles=0
		while (n<len(data)): #n is bits read so far
			t=struct.unpack('d',data[n:n+8])[0]
			Nsinks=struct.unpack('i',data[n+8:n+12])[0]
			
			#get sink masses
			m=0 #cumulative mass of all sinks
			for j in range(Nsinks):
				start = n + 8 + 4 + (j*bits_sink)
				m+=struct.unpack('d',data[start+8*9:start + 8*10])[0]

			#update arrays if Nsink changes
			if Nsinks != Nsinks_old:
				
				marker=0
				while marker == 0: #delete anything overlapping from last file 
					if len(T)>0:
						if t<T[-1]:
							print('overlap '+files[i]+' and '+files[i-1])
							T=np.delete(T,len(T)-1)
							N=np.delete(N,len(N)-1)
							M=np.delete(M,len(M)-1)
						else:
							marker=1
					else:
						marker=1
				#and update lists 
				N=np.append(N,Nsinks_old) #line to the right 
				N=np.append(N,Nsinks)     #and vertically up
				T=np.append(T,t)
				T=np.append(T,t)
				M=np.append(M,m)
				M=np.append(M,m)
			Nsinks_old = Nsinks

			#also update arrays every 10 cycles through the sink file (needed to track mass)
			if cycles>0:
				if Nsinks>0:
					if cycles/10==int(cycles/10):
						N=np.append(N,Nsinks)
						T=np.append(T,t)
						M=np.append(M,m)


			n  = n + 8 + 4 + Nsinks*bits_sink
			cycles+=1
				
	return T,N,M

def Nsink_plot(dirnames):
	fig,ax=plt.subplots(2,sharex=True)
	plt.subplots_adjust(hspace=0)
	colors='b','g','r','cyan','purple'
	labels=r'$\rho_{sink}$=10$^{-10}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-9}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-8}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-7}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-6}$gcm$^{-3}$'
	for i in range(len(dirnames)):
		t,N,M=Nsinks(dirnames[i])
		if i==0:
			t0=t[0]
		ax[0].plot((t-t0)*code_units.t_cu/(60*60*24*365),N,color=colors[i],label=labels[i])
		ax[1].plot((t-t0)*code_units.t_cu/(60*60*24*365),M*code_units.M_cu/ap.M_sun.cgs.value,color=colors[i],label=labels[i])
	ax[0].tick_params(axis="x", labelsize=10,direction="in",which='both')
	ax[0].tick_params(axis="y", labelsize=10,direction="in",which='both')
	ax[1].tick_params(axis="y", labelsize=10,direction="in",which='both')
	ax[1].tick_params(axis="x", labelsize=10,direction="in",which='both')
	ax[1].set_xlabel('t [yr]')
	ax[0].set_ylabel(r'N$_{\rm sink}$')
	ax[1].set_ylabel(r'M [M$_\odot$]')
	ax[1].legend(frameon=False,loc='lower right',fontsize=8)
	ax[0].set_ylim(0,N.max()+1)
	plt.show()
	return fig,ax

def Nsink_plot_big(dirnames):
	fig,ax=plt.subplots(ncols=3,nrows=2,sharex='col',sharey='row')
	plt.subplots_adjust(hspace=0,wspace=0,right=0.75)
	extensions='/1e8/sink_particle_info/','/1e9/sink_particle_info/','/1e10/sink_particle_info/','/1e11/sink_particle_info/','/1e12/sink_particle_info/'
	colors='b','g','r','cyan','purple'
	labels=r'$\rho_{sink}$=10$^{-10}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-9}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-8}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-7}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-6}$gcm$^{-3}$'
	for i in range(len(dirnames)):
		for j in range(len(extensions)):
			t,N,M=Nsinks(dirnames[i]+extensions[j])
			if j==0:
				t0=t[0]
			ax[0,i].plot((t-t0)*code_units.t_cu/(60*60*24*365),N,color=colors[j],label=labels[j])
			ax[1,i].plot((t-t0)*code_units.t_cu/(60*60*24*365),M*code_units.M_cu/ap.M_sun.cgs.value,color=colors[j],label=labels[j])
			ax[0,i].tick_params(axis="x", labelsize=10,direction="in",which='both')
			ax[0,i].tick_params(axis="y", labelsize=10,direction="in",which='both')
			ax[1,i].tick_params(axis="y", labelsize=10,direction="in",which='both')
			ax[1,i].tick_params(axis="x", labelsize=10,direction="in",which='both')
			ax[1,i].set_xlim(-10,1050)
			ax[1,i].set_xlabel('t [yr]',fontsize=10)
			ax[1,i].set_xticks([0,300,600,900])
			ax[0,i].set_title(('A','B','C')[i],fontsize=10)
		ax[0,0].set_ylabel(r'N$_{\rm sink}$',fontsize=10)
		ax[1,0].set_ylabel(r'M [M$_\odot$]',fontsize=10)	
		ax[0,-1].legend(frameon=False,loc=(1.05,0.3),fontsize=8)

def Nsink_MHD(dirnames):
	fig,ax=plt.subplots(ncols=3,nrows=2,sharey='row',sharex='col')
	plt.subplots_adjust(wspace=0,hspace=0)
	colors='aquamarine','royalblue','Crimson'
	labels='hydro',r'k$^{3/2}$','uniform'
	titles=r'10$^{-10}$gcm$^{-3}$',r'10$^{-9}$gcm$^{-3}$',r'10$^{-8}$gcm$^{-3}$'
	for i in range(len(dirnames)):
		t,N,M=Nsinks(dirnames[i]+'/sink_particle_info/')
		t0=t[0]
		ax[0,i].plot((t-t0)*code_units.t_cu/(60*60*24*365),N,color=colors[0],label=labels[0])
		ax[1,i].plot((t-t0)*code_units.t_cu/(60*60*24*365),M*code_units.M_cu/ap.M_sun.cgs.value,color=colors[0],label=labels[0])
		t,N,M=Nsinks(dirnames[i]+'MHD/sink_particle_info/')
		ax[0,i].plot((t-t0)*code_units.t_cu/(60*60*24*365),N,color=colors[1],label=labels[1])
		ax[1,i].plot((t-t0)*code_units.t_cu/(60*60*24*365),M*code_units.M_cu/ap.M_sun.cgs.value,color=colors[1],label=labels[1])
		t,N,M=Nsinks(dirnames[i]+'_uniform/sink_particle_info/')
		ax[0,i].plot((t-t0)*code_units.t_cu/(60*60*24*365),N,color=colors[2],label=labels[2])
		ax[1,i].plot((t-t0)*code_units.t_cu/(60*60*24*365),M*code_units.M_cu/ap.M_sun.cgs.value,color=colors[2],label=labels[2])
		ax[1,i].set_xlabel('t [yr]')
		ax[1,i].set_xlim(-10,2500)
		ax[0,i].set_title(titles[i])
	ax[-1,-1].legend(frameon=False,loc='lower right',fontsize=8)
	ax[0,0].set_ylabel(r'N$_{\rm sink}$')
	ax[1,0].set_ylabel(r'M [M$_\odot$]')
	return fig,ax
	
def IMF(dirname,t_goal,scale_ejects):
	'''if you don't care about ejections, set scale_ejects as 'no' '''
	files = allfiles(dirname)
	Nsinks_old =0
	M=np.array([])
	x=np.array([])
	y=np.array([])
	z=np.array([])
	mass_ejected=np.array([])
	t_old=-1
	for i in range(len(files)):
		if len(M)==0:
			with open(dirname+files[i], mode='rb') as file:
				data = file.read()
			bits_sink= 8*14 + 8 + 4*3 + 4
			n=0
			cycles=0
			while (n<len(data)): #n is bits read so far
				t=struct.unpack('d',data[n:n+8])[0]
				Nsinks=struct.unpack('i',data[n+8:n+12])[0]
				                 
                                               
				if t<t_old: #sometimes this happens when the sink files overlap
					t_old=t

				if abs(t-t_goal) > abs(t_old-t_goal): #reached the time goal
					for j in range(Nsinks):
						start = n + 8 + 4 + (j*bits_sink)
						Mcurrent=struct.unpack('d',data[start+8*9:start + 8*10])
						M=np.append(M,Mcurrent)
						
						if scale_ejects !='no':
							x=np.append(x,struct.unpack('d',data[start:start + 8]))
							y=np.append(y,struct.unpack('d',data[start+8:start + 8*2]))
							z=np.append(z,struct.unpack('d',data[start+8*2:start + 8*3]))
									
							
					n=len(data)

				t_old=t
				n  = n + 8 + 4 + Nsinks*bits_sink
				cycles+=1

	mid = np.where(M==M.max())
	if scale_ejects !='no':
		rs=np.sqrt((x-x[mid])**2+(y-y[mid])**2+(z-z[mid])**2)
		eject_mask=np.where(rs>scale_ejects)
		mass_ejected=M[eject_mask]
	
	return M,mass_ejected


def IMF_plot(dirnames,t_goal_yrs,scale_ejects):

	T,N,M=Nsinks(dirnames[0])
	t0=T[0]
	t_goal=(t_goal_yrs*60*60*24*365/code_units.t_cu)+t0
	m_min=1e12 * 4/3*np.pi*(1.71E-05)**3 *0.01
	rhos=r'$\rho_{sink}$=10$^{-10}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-9}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-8}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-7}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-6}$gcm$^{-3}$'
	colors='b','g','r','cyan','purple'
	minmass=np.array([1e8,1e9,1e10,1e11,1e12]) *code_units.rho_cu  * 4/3*np.pi * (np.array([0.001376823,0.0004563,0.000152667,5.04801E-05,1.71E-05])*code_units.d_cu)**3 /ap.M_sun.cgs.value
	fig,ax=plt.subplots(len(dirnames),sharex=True)
	plt.subplots_adjust(hspace=0)
	for i in range(len(dirnames)):
		M,Meject=IMF(dirnames[i],t_goal,scale_ejects)
		print(M)
		if i==0:
			m_max=M.max()*1.2
			bins = 10**np.linspace(np.log10(m_min),np.log10(m_max),100)
		ax[i].hist(M,bins,color=colors[i])
		if scale_ejects !='no':
			ax[i].hist(Meject,bins,color='k')
		#ax[i].axvline(x=minmass[i],ymin=0,ymax=1,color='k',linestyle='--')
		ax[i].set_xscale('log')
		N,M=np.histogram(M,bins)
		ax[i].set_ylim(0,N.max()+1)
		ax[i].set_yticks([N.max()])
		ax[i].text(1.22,0.5,rhos[i],ha='center', va='center', transform=ax[i].transAxes,fontsize=10)	
	ax[len(dirnames)-1].set_xlabel(r'M [M$_{\odot}$]')
	ax[int(len(dirnames)/2)].set_ylabel(r'N$_{\rm M}$             ',rotation=0)
	plt.subplots_adjust(left = 0.15,bottom = 0.17,right=0.7)
	return fig,ax

def IMF_plot_join(dirnames,t_goal_yrs,scale_ejects):
	fig,ax=plt.subplots(5,sharex=True)
	plt.subplots_adjust(hspace=0)
	colors='b','g','r','cyan','purple'
	rhos=r'$\rho_{sink}$=10$^{-10}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-9}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-8}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-7}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-6}$gcm$^{-3}$'
	extensions='1e8/sink_particle_info/','1e9/sink_particle_info/','1e10/sink_particle_info/','1e11/sink_particle_info/','1e12/sink_particle_info/'
	for j in range(5):
		Msinks=np.array([])
		ejects=np.array([])
		for i in range(len(dirnames)):
			T,N,M=Nsinks(dirnames[i]+extensions[0])
			t0=T[0]
			t_goal=(t_goal_yrs*60*60*24*365/code_units.t_cu)+t0
			Ms,Meject=IMF(dirnames[i]+extensions[j],t_goal,scale_ejects[i])
			print(Ms.min())
			Msinks=np.append(Msinks,Ms)
			if scale_ejects !='no':
				ejects=np.append(ejects,Meject)
		if j==0:
			m_max=M.max()*1.2
			m_min=1e12 * 4/3*np.pi*(1.71E-05)**3    *0.1
			bins = 10**np.linspace(np.log10(m_min),np.log10(m_max),100)
		ax[j].hist(Msinks,bins,color=colors[j])
		
		if scale_ejects !='no':
			ax[j].hist(ejects,bins,color='k')
		N,M=np.histogram(Msinks,bins)
		ax[j].set_ylim(0,N.max()+1)
		ax[j].set_yticks([N.max()])
		ax[j].text(1.22,0.5,rhos[j],ha='center', va='center', transform=ax[j].transAxes,fontsize=10)
		ax[j].set_xscale('log')
		ax[j].tick_params(axis="y", labelsize=10,direction="in",which='both')
		ax[j].tick_params(axis="x", labelsize=10,direction="in",which='both')
	plt.subplots_adjust(left = 0.15,bottom = 0.17,right=0.7)
	ax[-1].set_xlabel(r'M [M$_{\odot}$]')
	ax[2].set_ylabel(r'N$_{\rm M}$             ',rotation=0)

def IMF_plot_join_MHD(t_goal_yrs):

	fig,ax=plt.subplots(nrows=3,ncols=3,sharex='col')	
	plt.subplots_adjust(hspace=0,wspace=0)	
	dirname='/scratch/c.c1521474/resolution_test/'
	ext1='/MHD/','/MHD2/'
	ext2='/1e8','/1e9','/1e10'
	ext3='/sink_particle_info/','MHD/sink_particle_info/','_uniform/sink_particle_info/'	
	m_min=1e12 * 4/3*np.pi*(1.71E-05)**3    *0.1
	colors='b','g','r'

	for i in range(len(ext2)):
		M1=np.array([])
		M2=np.array([])
		M3=np.array([])
		for j in range(len(ext1)):
			t,N,M=Nsinks(dirname+ext1[j]+ext2[i]+'/sink_particle_info/')
			t0=t[0]
			t_goal=(np.array(t_goal_yrs)*60*60*24*365/code_units.t_cu)+t0
				
			Ms1,Meject1=IMF(dirname+ext1[j]+ext2[i]+ext3[0],t_goal[i],'no')
			Ms2,Meject2=IMF(dirname+ext1[j]+ext2[i]+ext3[1],t_goal[i],'no')
			Ms3,Meject3=IMF(dirname+ext1[j]+ext2[i]+ext3[2],t_goal[i],'no')

			M1=np.append(M1,Ms1)
			M2=np.append(M2,Ms2)
			M3=np.append(M3,Ms3)
		
		m_max=max(M1.max()*1.2,M2.max()*1.2,M3.max()*1.2)
		bins = 10**np.linspace(np.log10(m_min),np.log10(m_max),100)

		ax[0,i].hist(M1,bins,color=colors[i])
		ax[1,i].hist(M2,bins,color=colors[i])
		ax[2,i].hist(M3,bins,color=colors[i])
		ax[0,i].set_xscale('log')
		ax[1,i].set_xscale('log')
		ax[2,i].set_xscale('log')
			
def Nmerge(dirname):
	'''cycle starts wih double: time, int: Nsinks, then...
	doubles: 3pos, 3vel, 3accel, mass, FormationMass, FormationTime, AccretionRate, TimeOld =14
	long long ID, ints: HomeTask, Index, FormationOrder
	so sink mass starts at (8+4 + 8*9) bits, then every N(8*14 + 8 + 4*3 + 4) where N is the Number of sinks cycled through so far'''
	#first get all the sink_particle_info files into the right order
	files = allfiles(dirname)
	Nsinks_old =0
	T=np.array([])
	merge=np.array([])
	for i in range(len(files)):
		with open(dirname+files[i], mode='rb') as file:
			data = file.read()
		bits_sink= 8*14 + 8 + 4*3 + 4
		n=0
		cycles=0
		while (n<len(data)): #n is bits read so far
			t=struct.unpack('d',data[n:n+8])[0]
			Nsinks=struct.unpack('i',data[n+8:n+12])[0]
			if Nsinks < Nsinks_old:

				marker=0
				while marker == 0: #delete anything overlapping from last file
					if len(T)>0:
						if t<T[-1]:
							T=np.delete(T,len(T)-1)
							merge=np.delete(merge,len(merge)-1)
						else:
							marker=1
					else:	
						marker=1
				#and update lists
				if len(merge)==0:
					merge=np.append(merge,0) #line to the right
				else:
					merge=np.append(merge,merge[-1]+1)
				#merge=np.append(merge,merge[-1]+1)     #and vertically up
				#T=np.append(T,t)
				T=np.append(T,t)
			print(n,len(data)-( 8 + 4 + Nsinks*bits_sink))

			if n==(len(data)-( 8 + 4 + Nsinks*bits_sink)):
				if len(merge)==0:
					merge=np.append(merge,0)
				else:
					merge=np.append(merge,merge[-1])
				T=np.append(T,t)

			Nsinks_old=Nsinks 
			n  = n + 8 + 4 + Nsinks*bits_sink
			cycles+=1

	return T,merge	

def Nmerge_plot_join(dirnames):
	fig,ax=plt.subplots(1)
	plt.subplots_adjust(left=0.2,right=0.8,top=0.8,bottom=0.2)
	plt.subplots_adjust(hspace=0)
	extensions='1e8/sink_particle_info/','1e9/sink_particle_info/','1e10/sink_particle_info/','1e11/sink_particle_info/','1e12/sink_particle_info/'
	labels=r'$\rho_{sink}$=10$^{-10}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-9}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-8}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-7}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-6}$gcm$^{-3}$'
	colors='b','g','r','cyan','purple'
	for i in range(len(extensions)):
		
		MERGE=np.zeros(1000)
		for j in range(len(dirnames)):
			t,merge=Nmerge(dirnames[j]+extensions[i])
			t=t-t[0]
			if j==0:
				t_ = np.linspace(t[0],t[-2],1000)
			f=interp1d(t, merge,bounds_error=False,fill_value=merge.max())
			MERGE+=f(t_)
		ax.plot(t_ *code_units.t_cu /(60*60*24*365) , MERGE, c=colors[i],label=labels[i])
		ax.tick_params(axis="y", labelsize=10,direction="in",which='both')
		ax.tick_params(axis="x", labelsize=10,direction="in",which='both')
		ax.set_xlabel('t [yr]')
		ax.set_ylabel(r'N$_{\rm merge}$')
	ax.set_xlim(-10,400)
	ax.legend(frameon=False,fontsize=8)

def largest_sink(dirname):
	files = allfiles(dirname)
	M=np.array([])
	ACC=np.array([])
	T=np.array([])
	for i in range(len(files)):
		with open(dirname+files[i], mode='rb') as file:
			data = file.read()
		bits_sink= 8*14 + 8 + 4*3 + 4
		n=0
		cycles=0
		while (n<len(data)): #n is bits read so far
			#if cycles/interval==int(cycles/interval):
			t=struct.unpack('d',data[n:n+8])[0]
			T=np.append(T,t)

			Nsinks=struct.unpack('i',data[n+8:n+12])[0]
			mass=np.array([])			
			for j in range(Nsinks):
				start = n + 8 + 4 + (j*bits_sink)
				mass=np.append(mass,struct.unpack('d',data[start+8*9:start + 8*10]))
			
			M=np.append(M,mass.max())
			if cycles==0:
				ACC=np.append(ACC,0)
			else:
				acc=(mass.max()-mass_old) / ((t-t_old)*code_units.t_cu/(60*60*24*365)) #mass code units is already Msun
				ACC=np.append(ACC,acc)
			t_old=t
			mass_old=mass.max()
			n  = n + (8 + 4 + Nsinks*bits_sink) 
			cycles+=1
	return T,M,ACC

def largest_plot(dirnames,interval):
	colors='b','g','r','cyan','purple'
	labels=r'$\rho_{sink}$=10$^{-10}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-9}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-8}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-7}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-6}$gcm$^{-3}$'
	fig,ax=plt.subplots(2,sharex=True)
	plt.subplots_adjust(hspace=0)
	ax[1].set_ylabel(r'$\dot {\rm M}_{\rm largest}$ [M$_\odot$ yr$^{-1}$]',fontsize=10)
	ax[0].set_ylabel(r'M$_{\rm largest}$ [M$_\odot$]',fontsize=10)
	ax[1].set_xlabel('t [yr]',fontsize=10)
	for i in range(len(dirnames)):
		T,M,ACC=largest_sink(dirnames[i])
		if i==0:
			t0=T[0]
		T=(T-t0)*code_units.t_cu/(60*60*24*365)
		ax[0].semilogy(T,M,c=colors[i],label=labels[i])
		T=T[1::interval]
		ACC=ACC[1::interval] #acc units already cgs 
		ax[1].semilogy(T,ACC,c=colors[i])
	for i in range(2):
		ax[i].tick_params(axis="y", labelsize=10,direction="in",which='both')
		ax[i].tick_params(axis="x", labelsize=10,direction="in",which='both')
	ax[0].legend(frameon=False,loc='lower right',fontsize=8)
			
				
def largest_plot_big(dirnames):
	fig,ax=plt.subplots(ncols=3,nrows=2,sharex='col',sharey='row')
	plt.subplots_adjust(hspace=0,wspace=0,right=0.75)
	extensions='/1e8/sink_particle_info/','/1e9/sink_particle_info/','/1e10/sink_particle_info/','/1e11/sink_particle_info/','/1e12/sink_particle_info/'
	colors='b','g','r','cyan','purple'
	labels=r'$\rho_{sink}$=10$^{-10}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-9}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-8}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-7}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-6}$gcm$^{-3}$'
	for i in range(len(dirnames)):
		for j in range(len(extensions)):
			T,M,ACC=largest_sink(dirnames[i]+extensions[j])
			if j==0:
				t0=T[0]
			T=(T-t0)*code_units.t_cu/(60*60*24*365)
			ax[0,i].semilogy(T,M,c=colors[j],label=labels[j])
			T=T[1::50]		
			ACC=ACC[1::50] #acc units already cgs
			ax[1,i].semilogy(T,ACC,c=colors[j],label=labels[j])
			ax[1,i].set_xlabel('t [yr]',fontsize=10)
			ax[1,i].set_xlim(-50,1050)
		ax[0,i].set_title(('A','B','C')[i],fontsize=10)
		ax[1,i].set_xticks([0,300,600,900])
		for j in range(2):
			ax[j,i].tick_params(axis="y", labelsize=10,direction="in",which='both')
			ax[j,i].tick_params(axis="x", labelsize=10,direction="in",which='both')
	ax[0,-1].legend(frameon=False,loc=(1.05,0.3),fontsize=8)
	ax[1,0].set_ylabel(r'$\dot {\rm M}_{\rm largest}$ [M$_\odot$ yr$^{-1}$]',fontsize=10)
	ax[0,0].set_ylabel(r'M$_{\rm largest}$ [M$_\odot$]',fontsize=10)


'''estimate the number of ionising photons per second'''
def read_Hirano2014(txtfile):
	r=[]
	M=[]
	with open(txtfile) as f:
		for line in f.readlines():
			M.append(line.split()[0])
			r.append(line.split()[1])
	return np.asarray(M).astype(float) ,np.asarray(r).astype(float)

def photons(dirname):
	'''keeping it in base units this time'''
	#read in the mass and accretin rate onto the most massive sink
	t,M,ACC=largest_sink(dirname)

	#take R from mass-radius relationship (Hirano2014)
	M1,r1=read_Hirano2014('P2_curve.txt')
	f1=interp1d(M1,r1,fill_value='extrapolate')
	R=f1(M)

	#convert from solar values to base units
	ACC=ACC*ap.M_sun.value /(60*60*24*365)
	M=M*ap.M_sun.value
	R=R*ap.R_sun.value
	t=(t-t[0])*code_units.t_cu/(60*60*24*365) #years (only used in the graph)

	#accretion luminosity as function of time
	Lacc=M*ACC*ap.G.value/ R 
	
	#Temperature from Stefan-Boltzman
	T=(Lacc/(4*np.pi*R**2) /ap.sigma_sb.value)**(1/4)

	#Plank function
	v=10**np.linspace(np.log10(3.28e15),20,100) #Lynman limit and above 
	dv=np.zeros_like(v)
	dv[:-1]=v[1:]-v[:-1]
	
	Nion=np.array([])
	for i in range(len(T)):
		plank=2*ap.h.value*v**3/ap.c.value**2 *1/(np.exp(ap.h.value*v/(ap.k_B.value*T[i]))-1)
		N=sum(4*np.pi * plank*dv /(ap.h.value*v) * 4*np.pi*R[i]**2)
		Nion=np.append(Nion,N)

	return Nion,t,R,Lacc,T




'''first couple of functions are for getting sink info when they newly form'''


def readsink(filename,savefile,timemax):
	with open(filename, mode='rb') as file:
		data = file.read()
	#header: double Time, int Nsinks
	#doubles: 3pos, 3vel, 3accel, mass, FormationMass, FormationTime, AccretionRate, TimeOld =14 
	#long long ID, ints: HomeTask, Index, FormationOrder
	bits_sink= 8*14 + 8 + 4*3 + 4 #that 4 at the end... not sure what it is... think its a buffer... =32 every time
	
	#arrays for return 
	times=np.array([])
	dists=np.array([])

	X=[] 
	n=0
	i=0
	f=open(savefile,'a+')
	while (n<len(data)):
		t=struct.unpack('d',data[n:n+8])[0]
		Nsinks=struct.unpack('i',data[n+8:n+12])[0]
		if t<=timemax:
			if i>0:
				if Nsinks>Nsinks_old:
					print('new sink formed: '+str(Nsinks_old) + ' to ' + str(Nsinks))
					#temporary arrays 
					ids_old=np.array([])
					ids_new=np.array([])
					id_newsink=np.array([])
					arg_newsink=np.array([])
					x=np.array([])
					y=np.array([])
					z=np.array([])
					vx=np.array([])
					vy=np.array([])
					vz=np.array([])
					mass=np.array([])

					for j in range(Nsinks_old):
						#read in sink ids from the last iteration
						id_position = n_old + 8 + 4 + 8*14 #skip over the header and pos/vel/accel
						id_old=struct.unpack('q',data[id_position + (j*bits_sink) : id_position + (j*bits_sink) + 8])
						ids_old=np.append(ids_old,id_old)

					for j in range(Nsinks):
						#print('reading current sink data')
						start = n + 8 + 4 + (j*bits_sink)
						x=np.append(x,struct.unpack('d',data[start:start + 8]))
						y=np.append(y,struct.unpack('d',data[start+8:start + 8*2]))
						z=np.append(z,struct.unpack('d',data[start+8*2:start + 8*3]))
						vx=np.append(vx,struct.unpack('d',data[start+8*3:start + 8*4]))
						vy=np.append(vy,struct.unpack('d',data[start+8*4:start + 8*5]))
						vz=np.append(vz,struct.unpack('d',data[start+8*5:start + 8*6]))
						mass=np.append(mass,struct.unpack('d',data[start+8*9:start + 8*10]))
					
						id_position = n + 8 + 4 + 14*8 + (j*bits_sink)
						id_new=struct.unpack('q',data[id_position  : id_position  + 8])
						ids_new=np.append(ids_new,id_new)
						weird=struct.unpack('i',data[start  + bits_sink-4 : start  + bits_sink ])
						#print(weird)
						if id_new not in ids_old:
							#print('new sink id: '+str(id_new))
							id_newsink=np.append(id_newsink,id_new)
							arg_newsink=np.append(arg_newsink,j)
					#print('new ids' + str(ids_new))

					for j in range(len(arg_newsink)):
						arg=int(arg_newsink[j])
						dx=abs(x[arg] - x)
						dy=abs(y[arg] - y)
						dz=abs(z[arg] - z)
						r=np.sqrt(dx**2+dy**2+dz**2)
						rmin=np.sort(r)[1]
						arg_closest=int(np.where(r==np.sort(r)[1])[0])
						f.write(str(t) + ' ' + str(rmin) + ' ' + str(x[arg]) + ' ' + str(y[arg]) + ' ' + str(z[arg]) + ' ' +  
							str(vx[arg]) + ' ' + str(vy[arg]) + ' ' + str(vz[arg]) + ' ' + str(mass[arg]) + ' ' + 
							str(x[arg_closest]) + ' ' + str(y[arg_closest]) + ' ' + str(z[arg_closest]) + ' ' + 
							str(vx[arg_closest]) + ' ' + str(vy[arg_closest]) + ' ' + str(vz[arg_closest]) + ' ' +str(mass[arg_closest])
							+ '\n')
						times=np.append(times,t)
						dists=np.append(dists,rmin)
		i = i+1 
				
		Nsinks_old = Nsinks
		n_old = n 
		n  = n + 8 + 4 + Nsinks*bits_sink
	f.close()
	return times,dists


def read_sinkfile(txtfile):
	with open(txtfile) as f:
		t=np.array([])
		r=np.array([])
		x=np.array([])
		y=np.array([])
		z=np.array([])
		vx=np.array([])
		vy=np.array([])
		vz=np.array([])
		m=np.array([])
		x1=np.array([])
		y1=np.array([])
		z1=np.array([])
		vx1=np.array([])
		vy1=np.array([])
		vz1=np.array([])
		m1=np.array([])
		for line in f.readlines():
			t=np.append(t,line.split()[0])
			r=np.append(r,line.split()[1])
			x=np.append(x,line.split()[2])
			y=np.append(y,line.split()[3])
			z=np.append(z,line.split()[4])
			vx=np.append(vx,line.split()[5])
			vy=np.append(vy,line.split()[6])
			vz=np.append(vz,line.split()[7])
			m=np.append(m,line.split()[8])
			x1=np.append(x1,line.split()[9])
			y1=np.append(y1,line.split()[10])
			z1=np.append(z1,line.split()[11])
			vx1=np.append(vx1,line.split()[12])
			vy1=np.append(vy1,line.split()[13])
			vz1=np.append(vz1,line.split()[14])
			m1=np.append(m1,line.split()[15])
	return t.astype(float),r.astype(float),x.astype(float),y.astype(float),z.astype(float),vx.astype(float),vy.astype(float),vz.astype(float),m.astype(float),x1.astype(float),y1.astype(float),z1.astype(float),vx1.astype(float),vy1.astype(float),vz1.astype(float),m1.astype(float)


def hist_pannel(files):
	fig,axs=plt.subplots(len(files),sharex=True)
	plt.subplots_adjust(hspace=0)
	colors='b','g','r','cyan','purple'
	for i in range(len(files)):
		t,r,x,y,z,vx,vy,vz,m,x1,y1,z1,vx1,vy1,vz1,m1=read_sinkfile(files[-1-i])
		r=r*code_units.d_cu/ap.au.cgs.value
		if i==0:
			bins=10**np.linspace(np.log10(np.sort(r)[1]),np.log10(r.max()*2),100)
		axs[-1-i].hist(r,bins=bins,color=colors[-1-i])
	axs[0].set_xscale('log')
	axs[len(files)-1].set_xlabel(r'r$_{\rm min}$ [AU]')


def energy_check(files):
	for i in range(len(files)):
		t,r,x,y,z,vx,vy,vz,m,x1,y1,z1,vx1,vy1,vz1,m1=read_sinkfile(files[-1-i])
		V=np.sqrt((vx-vx1)**2+(vy-vy1)**2+(vz-vz1)**2) * code_units.v_cu
		M=m*code_units.M_cu
		M1=m1*code_units.M_cu
		R=r*code_units.d_cu
		Nbound=len(V[np.where(0.5*M*V**2 < ap.G.cgs.value*M*M1/R)])
		Nfree=len(V)-Nbound
		print('file '+str(-1-i)+': '+str(Nbound)+' bound, '+str(Nfree)+' free.')








'''next couple of funtions are for getting all sink info at all times'''

def get_all_info(dirname):
	sink_space = np.zeros((100, 3, 10000))
	sink_info_files = allfiles(dirname)
	for filename in sink_info_files: 
		with open(dirname+'/'+filename, mode='rb') as file:
			data = file.read()
		bits_sink= 8*14 + 8 + 4*3 + 4 #that 4 at the end... not sure what it is... think its a buffer... =32 every time
		X=[]
		n=0
		i=0
		mass_old=np.array([])
		ids_old=np.array([])
		while (n<len(data)):
			t=struct.unpack('d',data[n:n+8])[0]
		
			Nsinks=struct.unpack('i',data[n+8:n+12])[0]
	

			ids_keep=np.array([])
			mass_keep=np.array([])
			ids=np.array([])
			mass=np.array([])
			tag=0
			for j in range(Nsinks):
				start = n + 8 + 4 + (j*bits_sink)
				ids=np.append(ids,struct.unpack('q',data[start+8*14:start + 8*15]))
				mass=np.append(mass,struct.unpack('d',data[start+8*9:start + 8*10]))
			
			
				#begin storing the info... only read this if you really need to...
				#storage shape is 100,3,100
				#the first dimension allows for 100 possible sinks
				#second dimension gives mass, accretion rate and time
				#third dimension progresses the data in time (10000 time slots)

						
				if ids[j] not in sink_space[:,0,0]: #going to put the info into the sink space array, each sink gets its own 2D array shape (3,10000)
					tag+=1
					I=0 #I is the arg of possible sink spaces in sink_space
					marker=0
					while marker==0:					
						if sink_space[I,0,0]==0: #found an unoccupied 2D array
							sink_space[I,0,0]=ids[j] #left with 2D array, first row contains id and number of entries
							sink_space[I,1,0]=2
						
							sink_space[I,0,1]=mass[j] #second row onwards contains the entries of mass + accretion rate 
							sink_space[I,1,1]=0 #initial accretion = 0 bc no previous data 
							sink_space[I,2,1]=t 
							marker=1 #let the loop stop looking fo a space 
						I+=1

				if ids[j] in sink_space[:,0,0]:
					tag+=1
					arg = np.where(sink_space[:,0,0]==ids[j])[0][0]
					entry = int(sink_space[arg,1,0])
				
					sink_space[arg,0,entry] = mass[j]
				
					sink_space[arg,1,entry] = mass[j] - (sink_space[arg,0,entry-1]) / (sink_space[arg,2,entry-1]-t)
					sink_space[arg,2,entry]=t
					sink_space[arg,1,0]+=1 #update how long the list is 
				
				#ok the storing is done, I hope
				
				ids_old = ids
				mass_old = mass
			
			n  = n + 8 + 4 + Nsinks*bits_sink
			i+=1

	#time to read sink_space and plot it
	fig,ax=plt.subplots(1)
#	fig1,ax1=plt.subplots(1)
	for i in range(len(sink_space[:,0,0])-1):
		if sink_space[i,0,0]>0:
			ids=sink_space[i,0,0]
			number=int(sink_space[i,1,0])
			mass=sink_space[i,0,1:number+1] * code_units.M_cu / ap.M_sun.cgs.value
			acc=sink_space[i,1,1:number+1] * code_units.M_cu / ap.M_sun.cgs.value / (code_units.t_cu / (60*60*24*365))
			t=sink_space[i,2,1:number+1] * (code_units.t_cu / (60*60*24*365))
			mask=np.where(acc>0)
			ax.loglog(t,mass,label=int(ids))
	ax.legend(fontsize=6,loc=(1,-0.2))
	plt.subplots_adjust(right=0.75)
	
	return sink_space

def triple_mass_plot(dirnames):
	fig,ax=plt.subplots(ncols=3,sharey=True)
	plt.subplots_adjust(wspace=0)
	
	for j in range(3):
		ax[j].set_title(('A','B','C')[j],fontsize=10)
		ax[j].set_xlabel('t [yr]')
		ax[j].tick_params(axis="x", labelsize=10,direction="in",which='both')
		ax[j].tick_params(axis="y", labelsize=10,direction="in",which='both')
	ax[0].set_ylabel('M [M$_\odot$]')
	for j in range(len(dirnames)):
		sink_space=get_all_info(dirnames[j])
		for i in range(len(sink_space[:,0,0])-1):
			if sink_space[i,0,0]>0:
				ids=sink_space[i,0,0]
				number=int(sink_space[i,1,0])
				mass=sink_space[i,0,1:number+1] * code_units.M_cu / ap.M_sun.cgs.value
				acc=sink_space[i,1,1:number+1] * code_units.M_cu / ap.M_sun.cgs.value / (code_units.t_cu / (60*60*24*365))
				t=sink_space[i,2,1:number+1] * (code_units.t_cu / (60*60*24*365))
				t=t[:-1]#last t is messed up, not sure why
				if i==0:
					t0=t[0]
				mass=mass[:-1]
				mask=np.where(acc>0)
				ax[j].semilogy(t-t0,mass,c='k')
				#ax[j].set_xlim((t-t[0]).min(),(t-t[0]).max())
