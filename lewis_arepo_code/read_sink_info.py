import numpy as np 
import matplotlib.pyplot as plt 
import struct
import code_units
import astropy.constants as ap 
import os

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
	files = allfiles(dirname)
	Nsinks_old =0 
	T=np.array([])
	N=np.array([])
	for i in range(len(files)):
		
		with open(dirname+files[i], mode='rb') as file:
			data = file.read()
		bits_sink= 8*14 + 8 + 4*3 + 4
		n=0
		while (n<len(data)):
			t=struct.unpack('d',data[n:n+8])[0]
			Nsinks=struct.unpack('i',data[n+8:n+12])[0]
			if Nsinks != Nsinks_old:
				if len(T>0):
					while t<T[-1]:
						T=np.delete(T,len(T)-1)
						T=np.delete(T,len(T)-1)
						N=np.delete(N,len(N)-1)
						N=np.delete(N,len(N)-1)
				N=np.append(N,Nsinks_old) #line to the right 
				N=np.append(N,Nsinks)     #and vertically up
				T=np.append(T,t)
				T=np.append(T,t)
				Nsinks_old = Nsinks 
			n  = n + 8 + 4 + Nsinks*bits_sink
	return T,N

def Nsink_plot(dirnames):
	colors='b','g','r','cyan','purple'
	for i in range(len(dirnames)):
		t,N=Nsinks(dirnames[i])
		plt.plot(t,N,color=colors[i])
	plt.show()


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

def get_all_info(sink_info_files,savefile):
	sink_space = np.zeros((100, 3, 10000))
	for filename in sink_info_files: 
		with open(filename, mode='rb') as file:
			data = file.read()
		bits_sink= 8*14 + 8 + 4*3 + 4 #that 4 at the end... not sure what it is... think its a buffer... =32 every time
		X=[]
		n=0
		i=0
		mass_old=np.array([])
		ids_old=np.array([])
		entries=0
		while (n<len(data)):
			t=struct.unpack('d',data[n:n+8])[0]
			Nsinks=struct.unpack('i',data[n+8:n+12])[0]
			entries+=1
			if entries == 1:
				t0=t

			ids_keep=np.array([])
			mass_keep=np.array([])
			ids=np.array([])
			mass=np.array([])
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
					arg = np.where(sink_space[:,0,0]==ids[j])[0][0]
					entry = int(sink_space[arg,1,0])
				
					sink_space[arg,0,entry] = mass[j]
				
					sink_space[arg,1,entry] = (mass[j] - sink_space[arg,0,entry-1]) / (t-sink_space[arg,2,entry-1])
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
	starttime=sink_space[0,2,1]
	for i in range(len(sink_space[:,0,0])-1):
		if sink_space[i,0,0]>0:
			ids=sink_space[i,0,0]
			number=int(sink_space[i,1,0])
			mass=sink_space[i,0,1:number+1] * code_units.M_cu / ap.M_sun.cgs.value
			acc=sink_space[i,1,1:number+1] * code_units.M_cu / ap.M_sun.cgs.value / (code_units.t_cu / (60*60*24*365))
			t=(sink_space[i,2,1:number+1]-starttime) * (code_units.t_cu / (60*60*24*365)) 
			mask=np.where(acc>0)
			ax.loglog(mass[np.where(acc>0)],acc[np.where(acc>0)],label=int(ids))
#	ax.legend(fontsize=6,loc=(1,-0.2))
#	plt.subplots_adjust(right=0.75)
	ax.tick_params(axis="x", labelsize=10,direction="in",which='both')
	ax.tick_params(axis="y", labelsize=10,direction="in",which='both')

