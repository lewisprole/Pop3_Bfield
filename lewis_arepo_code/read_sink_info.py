import numpy as np 
import matplotlib.pyplot as plt 
import struct
import code_units
import astropy.constants as ap 


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

					for j in range(Nsinks_old):
						#print('checking old ids')
						#read in sink ids from the last iteration
						id_position = n_old + 8 + 4 + 8*14 #skip over the header and pos/vel/accel
						id_old=struct.unpack('q',data[id_position + (j*bits_sink) : id_position + (j*bits_sink) + 8])
						ids_old=np.append(ids_old,id_old)
					#print('old ids' + str(ids_old))

					for j in range(Nsinks):
						#print('reading current sink data')
						start = n + 8 + 4 + (j*bits_sink)
						x=np.append(x,struct.unpack('d',data[start:start + 8]))
						y=np.append(y,struct.unpack('d',data[start+8:start + 8*2]))
						z=np.append(z,struct.unpack('d',data[start+8*2:start + 8*3]))
					
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
						f.write(str(t) + ' ' + str(rmin) + ' ' + str(x[arg]) + ' ' + str(y[arg]) + ' ' + str(z[arg]) + '\n')
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
		t=[]
		r=[]
		x=[]
		y=[]
		z=[]
		for line in f.readlines():
			t.append(line.split()[0])
			r.append(line.split()[1])
			x.append(line.split()[2])
			y.append(line.split()[3])
			z.append(line.split()[4])
	return t,r,x,y,z


def hist_pannel(files):
	fig,axs=plt.subplots(len(files),sharex=True)
	plt.subplots_adjust(hspace=0)
	colors='b','g','r','cyan','purple'
	for i in range(len(files)):
		t,r,x,y,z=read_sinkfile(files[-1-i])
		r=np.asarray(r).astype(float)*code_units.d_cu/ap.au.cgs.value
		if i==0:
			bins=10**np.linspace(np.log10(np.sort(r)[1]),np.log10(r.max()*2),100)
		axs[-1-i].hist(r,bins=bins,color=colors[-1-i])
	axs[0].set_xscale('log')
	axs[len(files)-1].set_xlabel(r'r$_{\rm min}$ [AU]')
