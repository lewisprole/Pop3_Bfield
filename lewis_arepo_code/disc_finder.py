import numpy as np 
import matplotlib.pyplot as plt 
import sink_mask
import code_units 
import arepo_utils
import astropy.constants as ap
from matplotlib.colors import LogNorm
import read_sink_info  


def another_try(a,rho_sink,zoomzone):
	'''this time stay in spherical coordinates, care about density drop AND large vrot AND bound'''
	sinks=[]
	for i in range(a.npart[-1]):
		sinks.append([])

	bound_to=np.ones(len(a.x)) * -1 
	mask=sink_mask.sinkmask(a,zoomzone)[0]
	M=a.mass * code_units.M_cu
	Msink=a.sinkmass *code_units.M_cu
	for i in range(len(a.x[mask])):
		r=np.sqrt((a.x[mask[i]]-a.sinkx)**2+(a.y[mask[i]]-a.sinky)**2+(a.z[mask[i]]-a.sinkz)**2) *code_units.d_cu
		forces=ap.G.cgs.value * M[mask[i]] * Msink  /r
		vmag=np.sqrt(a.vx[mask[i]]**2+a.vy[mask[i]]**2+a.vz[mask[i]]**2) *code_units.v_cu
		if 0.5 * M[mask[i]] *vmag**2 < forces.max():
			bound_to[mask[i]]=np.where(forces==forces.max())[0]

	for i in range(int(a.npart[-1])):
		mask_sink=np.where(bound_to==i)[0] #the cells bound to this sink 
		r=np.sqrt((a.x[mask_sink]-a.sinkx[i])**2+(a.y[mask_sink]-a.sinky[i])**2+(a.z[mask_sink]-a.sinkz[i])**2) *code_units.d_cu
		args=np.argsort(r)
		rho_old = rho_sink
		for J in range(len(args)):
			j=mask_sink[args[J]]
			if a.rho[j]>rho_old/10:
				rvec=np.array([ a.x[j]-a.sinkx[i], a.y[j]-a.sinky[i], a.z[j]-a.sinkz[i] ]) *code_units.d_cu
				vvec=np.array([ a.vx[j], a.vy[j], a.vz[j] ]) * code_units.v_cu
				vmag=np.sqrt(a.vx[j]**2+a.vy[j]**2+a.vz[j]**2) * code_units.v_cu

				cross1 = vvec[1] * rvec[2] - vvec[2] * rvec[1]
				cross2 = -vvec[0] * rvec[2] + vvec[2] * rvec[0]
				cross3 = vvec[0] * rvec[1] - vvec[1] * rvec[0]
				vrot = 1/r[args[J]] * np.array([cross1,cross2,cross3]) #velocity perpendicular to r (to sink)
				vrotmag=np.sqrt(vrot[0]**2+vrot[1]**2+vrot[2]**2)
				if vrotmag> 0.5*vmag:
					sinks[i].append(j)
					rho_old = a.rho[j]
				
				
				

'''THIS IS THE ONE'''



def fifty_close(a,sinkx,sinky,sinkz,accradius):
	'''look at 50 cells closest to sink
	assume that cells in disk have higher ratio of rotational velocity to total velocity 
	take the unit vector of the rotational velocity to be a rod through the disk
	return the rotational velocity unit vector
	'''
	r=np.sqrt((a.x-sinkx)**2+(a.y-sinky)**2+(a.z-sinkz)**2) *code_units.d_cu
	mask=np.where(r<accradius*2 *code_units.d_cu)[0] #new
	if len(mask)<50:
		r_args=np.argsort(r)
		mask=r_args[:50] #50 close lol
	args=np.argsort(r[mask]) #new
	#mask=args[0:50]
	RATIO=np.array([])
	ratiomax=0
	X=np.array([])
	Y=np.array([])
	Z=np.array([])
	#for i in range(len(mask)):
		#I=mask[i]
	for i in range(len(args)):
		I=mask[args[i]]
		rvec=np.array([ a.x[I]-sinkx, a.y[I]-sinky, a.z[I]-sinkz]) *code_units.d_cu
		vvec=np.array([ a.vx[I], a.vy[I], a.vz[I] ]) * code_units.v_cu
		cross1 = vvec[1] * rvec[2] - vvec[2] * rvec[1]
		cross2 = -vvec[0] * rvec[2] + vvec[2] * rvec[0]
		cross3 = vvec[0] * rvec[1] - vvec[1] * rvec[0]
		vrot = 1/r[I] * np.array([cross1,cross2,cross3]) #velocity perpendicular to r (to sink) 

		vrotmag=np.sqrt(vrot[0]**2 + vrot[1]**2 + vrot[2]**2)
		vrotunit=vrot/vrotmag #points perpendicular to the disk
		vmag=np.sqrt(a.vx[I]**2+ a.vy[I]**2+ a.vz[I]**2)* code_units.v_cu
		ratio = vrotmag/vmag #ratio of rotational velocity to total velocity
		if ratio>0.6:
			X=np.append(X,vrotunit[0])
			Y=np.append(Y,vrotunit[1])
			Z=np.append(Z,vrotunit[2])
	x=np.median(X)
	y=np.median(Y)
	z=np.median(Z)
	vector = np.array([x,y,z])
	mag=np.sqrt(x**2+y**2+z**2) #probably not quite a unit vector yet because of the averaging
	vector=vector/mag

	return vector 

def combine_binaries(a,merge_length,accradius):
	#mergers 
	merge=np.linspace(0,a.npart[-1]-1,a.npart[-1]).astype(int) #args of array correspond to the sinki
	for i in range(a.npart[-1]):
		if merge[i]==i:
			r=np.sqrt((a.sinkx[i]-a.sinkx)**2+(a.sinky[i]-a.sinky)**2+(a.sinkz[i]-a.sinkz)**2)
			to_merge=np.where(r<merge_length)[0]
			if len(to_merge)>0:
				for j in range(len(to_merge)): #any sinks previously merging with to_merge[j] are now going to i 
					merge[np.where(merge==to_merge[j])]=i
	logged=np.array([])
	sinkx=np.array([])
	sinky=np.array([])
	sinkz=np.array([])
	sinkmass=np.array([])
	sinkvx=np.array([])
	sinkvy=np.array([])
	sinkvz=np.array([])
	ids=np.array([])
	vectors=[]
	for i in range(len(merge)):
		if merge[i] not in logged:
			combine=np.where(merge==merge[i])[0]
			ids=np.append(ids,a.sinkid[combine][np.where(a.sinkmass[combine]==a.sinkmass[combine].max())])
			x=sum(a.sinkmass[combine]*a.sinkx[combine])/sum(a.sinkmass[combine])
			y=sum(a.sinkmass[combine]*a.sinky[combine])/sum(a.sinkmass[combine])
			z=sum(a.sinkmass[combine]*a.sinkz[combine])/sum(a.sinkmass[combine])
			mass=sum(a.sinkmass[combine])
			vx=sum(a.sinkvx[combine])
			vy=sum(a.sinkvy[combine])
			vz=sum(a.sinkvz[combine])
	
			sinkx=np.append(sinkx,x)
			sinky=np.append(sinky,y)
			sinkz=np.append(sinkz,z)
			sinkmass=np.append(sinkmass,mass)
			sinkvx=np.append(sinkvx,vx)
			sinkvy=np.append(sinkvy,vy)
			sinkvz=np.append(sinkvz,vz)
			logged=np.append(logged,merge[i])

			#sort the vectors 
			vecx=np.array([])
			vecy=np.array([])
			vecz=np.array([])
			for j in range(len(combine)):
				vector=fifty_close(a,a.sinkx[combine][j],a.sinky[combine][j],a.sinkz[combine][j],accradius)
				vecx=np.append(vecx,vector[0])
				vecy=np.append(vecy,vector[1])
				vecz=np.append(vecz,vector[2])
			vector=np.array([np.mean(vecx),np.mean(vecy),np.mean(vecz)])
			vectors.append(vector)
				

	return sinkx,sinky,sinkz,sinkmass,sinkvx,sinkvy,sinkvz,ids,vectors


'''MASS WIEGHT EVERYTHING - sound speed, angular velocity for T
make combines IMF and IMF just for ejected
switch off hydro timebins -tree based timesteps off, force equal timesteps on 
read domain options in param files (arepo pdf)
make MHD IMFs (combined)'''	
def radial_Q(vel_vec,r_vec,temp,mass,Msink):

	r=np.sqrt((r_vec[0])**2+(r_vec[1])**2+(r_vec[2])**2)
	

	rs=np.linspace(r.min(),r.max(),50)
	Q=np.array([])
	for i in range(len(rs)-1):	
		mask=np.where((r>rs[i]) & (r<rs[i+1]))
		if len(mask[0])>0:
			dr=rs[i+1]-rs[i]
			dA=2*np.pi*rs[i] *dr
			T_massweighted=np.sum(temp[mask]*mass[mask])/np.sum(mass[mask])
			c_s=np.mean(np.sqrt(ap.k_B.cgs.value * T_massweighted/ap.m_p.cgs.value))
			surface=np.sum(mass[mask]) /dA
			
			vvec=vel_vec[0][mask], vel_vec[1][mask], vel_vec[2][mask] #strip down the vectors to within anulus
			rvec=r_vec[0][mask],r_vec[1][mask],r_vec[2][mask]

			crossi= vvec[1]*rvec[2] - vvec[2]*rvec[1] #rotational velocity calculation 
			crossj= -vvec[0]*rvec[2] + vvec[2]*rvec[0]
			crossz=vvec[0]*rvec[1] - vvec[1]*rvec[0]
			v_cross_r=np.sqrt(crossi**2+crossj**2+crossz**2)
			vr=v_cross_r / r[mask]
			vr_massweighted=np.sum(vr*mass[mask])/np.sum(mass[mask])

			circumference=2*np.pi*rs[i] #keplarian period 
			period=circumference/vr_massweighted
			#period=2*np.pi *np.sqrt(rs[i]**3 /(ap.G.cgs.value*Msink))
			frequency=2*np.pi/period 

			Q=np.append(Q,c_s*frequency/(np.pi*ap.G.cgs.value*surface))
		else:
			Q=np.append(Q,np.nan)
	return Q,rs

		
		
		
def get_disc(a,accradius,zoomzone,maxsize,merge_length):
	'''First merge the sinks that are part of the same disc, set by merge_length.
	Within radius maxsize of each sink, energy check to see what is gravitationally bound.
	Find vector through the disc using the rotational velocity vector.
	Filter out particles above 45deg from plane of disc.
	Filer out particles below 1/10th of running mean density of disc.
	Calculate radial profile of Toomr parameter Q.
	Return minimum Q, time, and lead sink ID of the disc (to track the disc between snapshots)'''
	
	Qmins=np.array([])
	sinkx,sinky,sinkz,sinkmass,sinkvx,sinkvy,sinkvz,ids,vectors=combine_binaries(a,merge_length,accradius) #merging 
	#print(1)
	sinks=[] #create array for each of the sinks, holds index of cells in their disc 
	for i in range(len(sinkmass)):
		sinks.append([])

	bound_to=np.ones(len(a.x)) * -1 #array of -1's, to be set to the sinki of the sink they are bound to
	mask=sink_mask.sinkmask(a,zoomzone)[0] #only do this in the very central region of the simulation box 
	M=a.mass * code_units.M_cu
	Msink=sinkmass *code_units.M_cu
	#print(2)
	for i in range(len(a.x[mask])): #loop through cells to find who they are bound to 
		r=np.sqrt((a.x[mask[i]]-sinkx)**2+(a.y[mask[i]]-sinky)**2+(a.z[mask[i]]-sinkz)**2) *code_units.d_cu #to all sinks
		forces=ap.G.cgs.value * M[mask[i]] * Msink  /r #grav force to all sinks
		vmag=np.sqrt(a.vx[mask[i]]**2+a.vy[mask[i]]**2+a.vz[mask[i]]**2) *code_units.v_cu

		found_sink=0
		j=0
		order=np.argsort(forces)
		while found_sink==0: #only bother assigning it to the sink if it's within the max disc radius
			sink_index=-1
			if r[order[-1-j]] <= maxsize*code_units.d_cu:
				sink_index=order[-1-j]
				found_sink=1
			j+=1
			if j>len(order)-1:
				found_sink=1
			
		if sink_index>-1:
			if 0.5 * M[mask[i]] *vmag**2 < 1.5*forces[sink_index]: #if it's bound,check it's in circular motion
				#rvec=np.array([(a.x[mask[i]]-a.sinkx[sink_index]),(a.y[mask[i]]-a.sinky[sink_index]),(a.z[mask[i]]-a.sinkz[sink_index])]) *code_units.d_cu
				#vvec=np.array([a.x[mask[i]],a.y[mask[i]],a.z[mask[i]]]) *code_units.v_cu
				#cross1 = vvec[1] * rvec[2] - vvec[2] * rvec[1]
				#cross2 = -vvec[0] * rvec[2] + vvec[2] * rvec[0]
				#cross3 = vvec[0] * rvec[1] - vvec[1] * rvec[0]
				#vrot = 1/r[sink_index] * np.sqrt(cross1**2+cross2**2+cross3**2) #velocity perpendicular to r (to sink)
				#if vrot>0.5*vmag: #more than half the velocity is circular	
				bound_to[mask[i]]=sink_index #index of the sink they are bound to
	#print(3)
	#fig1,ax1=plt.subplots(1)
	#for i in range(len(sinkmass)):
	#	bound=np.where(bound_to==i)[0]
	#	ax1.scatter(a.x[bound],a.y[bound],s=1)
	#print(4)

	#fig2,ax2=plt.subplots(1)
	#X,Y,Z=np.histogram2d( a.y[mask],a.x[mask],bins=(500,500))
	#Size,Y,Z=np.histogram2d( a.y[mask],a.x[mask],weights=(a.rho[mask])*code_units.rho_cu,bins=(500,500))
	#ax2.imshow(np.log10(Size/X),cmap='magma',aspect='auto',extent=[Z[0],Z[-1],Y[-1],Y[0]])

	#print(5)
	#arepo_utils.arepoimage(a.x[mask],a.y[mask],a.rho[mask])
	#fig3,ax3=plt.subplots(1)
	for i in range(len(sinkmass)): #loop through each sink
		print('sinking')
		vector=vectors[i]#fifty_close(a,sinkx[i],sinky[i],sinkz[i],accradius) #find vector through the disc 
		bound=np.where(bound_to==i)[0]
		dx=(a.x[bound]-sinkx[i]) *code_units.d_cu
		dy=(a.y[bound]-sinky[i]) *code_units.d_cu
		dz=(a.z[bound]-sinkz[i]) *code_units.d_cu
		r=np.sqrt(dx**2+dy**2+dz**2)
		vecmag=np.sqrt(vector[0]**2 + vector[1]**2 +vector[2]**2)
		costheta = (vector[0]*dx +  vector[1]*dy +  vector[2]*dz) / r / vecmag
		angle=np.where(np.arccos(costheta)>0.5*np.pi/4)[0]
		
		args=np.argsort(r[angle])
		marker=0
		j=0
		#find disc cells 
		while marker==0:
			if j<len(args):
				if j==0:
					rho_av=a.rho[bound][angle][args[j]]
					sinks[i].append(bound[angle][args[j]])
				
				else:
					if r[angle][args[j]]<maxsize*code_units.d_cu:
						if a.rho[bound][angle][args[j]]>rho_av/20:
							sinks[i].append(bound[angle][args[j]])
							rho_av=np.mean(a.rho[sinks[i]])
				j+=1
			else:
				marker=1
		print('toomring')
		#calculate Toomr parameter 
		if len(sinks[i])>0:
			rvec=np.array([a.x[sinks[i]]-sinkx[i] , a.y[sinks[i]]-sinky[i], a.z[sinks[i]]-sinkz[i] ])*code_units.d_cu
			vvec=np.array([a.vx[sinks[i]],a.vy[sinks[i]],a.vz[sinks[i]]]) * code_units.v_cu
			Qs,R=radial_Q(vvec,rvec,a.temp[sinks[i]],a.mass[sinks[i]]*code_units.M_cu,sinkmass[i]*code_units.M_cu)
			#R=np.sqrt((a.x[sinks[i]]-a.sinkx[i])**2+(a.y[sinks[i]]-a.sinky[i])**2+(a.z[sinks[i]]-a.sinkz[i])**2)
			#ax3.semilogy(R[:-1]/ap.au.cgs.value,Qs,'.')
			#ax2.scatter(a.x[sinks[i]],a.y[sinks[i]],s=0.1)
			if len(Qs[~np.isnan(Qs)])>0:
				Qmins=np.append(Qmins,Qs[~np.isnan(Qs)].min())
			else:	
				Qmins=np.append(Qmins,np.nan)
		else:
			Qmins=np.append(Qmins,np.nan)
		
	print(6)
	#ax3.axhline(y=1,c='k'),plt.xlabel('R [AU]'),plt.ylabel('Q')
	#ax2.plot(a.sinkx,a.sinky,'x',c='k')
	#ax2.plot(sinkx,sinky,'o',c='k')

	return Qmins,ids,a.time



def Q_time(dirname,start,end,interval,accradius,zoomzone,maxsize,merge_length):
	Qs=[]
	T=[]
	IDS=np.array([])
	#Nsinks=np.array([])
	#Nsinkst=np.array([])
	for i in range(int((end-start)/interval)):
		n=str(start+i*interval)
		if (start+i*interval)<100:
			n='0'+str(start+i*interval)
			if (start+i*interval)<10:
				n='00'+str(start+i*interval)
		a=arepo_utils.aread(dirname+'/snapshot_'+n)
		#if i==0:
			#t0=a.time
		#	Nsink=a.npart[-1]
		Qmins,ids,t=get_disc(a,accradius,zoomzone,maxsize,merge_length)
		for j in range(len(Qmins)):
			if ids[j] not in IDS:
				IDS=np.append(IDS,ids[j])
				Qs.append([])
				T.append([])
				index=len(IDS)-1
				Qs[index].append(Qmins[j])
				T[index].append(t)
			else:
				index=np.where(IDS==ids[j])[0][0]
				Qs[index].append(Qmins[j])
				T[index].append(t)
		#if a.npart[-1]>Nsink:
		#	Nsinks=np.append(Nsinks,a.npart[-1])
		#	Nsinkst=np.append(Nsinkst,t)
		#Nsink=a.npart[-1]
	fig,ax=plt.subplots(1)
	sinktime,N,M=read_sink_info.Nsinks(dirname+'/sink_particle_info/')
	#plt.figure()
	Nsink=N[0]
	t0=sinktime[0]
	for i in range(len(sinktime)):
		if N[i]>Nsink:
			plt.axvline(x=(sinktime[i]-t0)*code_units.t_cu/(60*60*24*365),c='r')
		Nsink=N[i]
	for i in range(len(Qs)):
 		ax.semilogy((T[i]-t0)*code_units.t_cu/(60*60*24*365),Qs[i],c='k')
	
		
	return Qs,T,IDS
	
def Qjoin(dirnames):
	fig,ax=plt.subplots(nrows=len(dirnames),sharex=True)
	plt.subplots_adjust(hspace=0)
	rhos=r'$\rho_{\rm sink}$=10$^{-10}$gcm$^{-3}$',r'$\rho_{\rm sink}$=10$^{-9}$gcm$^{-3}$',r'$\rho_{\rm sink}$=10$^{-8}$gcm$^{-3}$',r'$\rho_{\rm sink}$=10$^{-7}$gcm$^{-3}$',r'$\rho_{\rm sink}$=10$^{-6}$gcm$^{-3}$'
	for j in range(len(dirnames)):
		sinktime,N,M=read_sink_info.Nsinks(dirnames[j]+'/sink_particle_info/')
		Nsink=N[0]
		t0=sinktime[0]
		for i in range(len(sinktime)):
			if N[i]>Nsink:
				ax[j].axvline(x=(sinktime[i]-t0)*code_units.t_cu/(60*60*24*365),c='r')
			Nsink=N[i]
		if j==0:
			Qs,T,IDS=Q_time('/scratch/c.c1521474/resolution_test/merge/1e8_redo/',43,143,5,0.001376823,0.01,0.0025,0.0001)
			for k in range(len(Qs)):
				ax[j].semilogy((T[k]-t0)*code_units.t_cu/(60*60*24*365),Qs[k],c='k')
		if j==1:
			Qs,T,IDS=Q_time('/scratch/c.c1521474/resolution_test/merge/1e9/',36,146,5,0.001376823,0.01,0.0025,0.0001)
			for k in range(len(Qs)):
				ax[j].semilogy((T[k]-t0)*code_units.t_cu/(60*60*24*365),Qs[k],c='k')
		if j==2:
			Qs,T,IDS=Q_time('/scratch/c.c1521474/resolution_test/merge/1e10/',34,144,5,0.001376823,0.01,0.0025,0.0001)
			for k in range(len(Qs)):
				ax[j].semilogy((T[k]-t0)*code_units.t_cu/(60*60*24*365),Qs[k],c='k')
		if j==3:
			Qs,T,IDS=Q_time('/scratch/c.c1521474/resolution_test/merge/1e11/',38,148,5,0.001376823,0.01,0.0025,0.0005)
			for k in range(len(Qs)):
				ax[j].semilogy((T[k]-t0)*code_units.t_cu/(60*60*24*365),Qs[k],c='k')
		if j==4:
			Qs,T,IDS=Q_time('/scratch/c.c1521474/resolution_test/merge/1e12_redo/',45,185,10,1.71E-05,0.005,0.0002,0.0001)
			for k in range(len(Qs)):
				ax[j].semilogy((T[k]-t0)*code_units.t_cu/(60*60*24*365),Qs[k],c='k')
			Qs,T,IDS=Q_time('/scratch/c.c1521474/resolution_test/merge/1e12_redo/',195,355,10,1.71E-05,0.01,0.005,0.001)
			for k in range(len(Qs)):
				ax[j].semilogy((T[k]-t0)*code_units.t_cu/(60*60*24*365),Qs[k],c='k')
			Qs,T,IDS=Q_time('/scratch/c.c1521474/resolution_test/merge/1e12_redo/',355,545,10,1.71E-05,0.013,0.005,0.001)	
			for k in range(len(Qs)):
				ax[j].semilogy((T[k]-t0)*code_units.t_cu/(60*60*24*365),Qs[k],c='k')

		ax[j].text(1.22,0.5,rhos[j],ha='center', va='center', transform=ax[j].transAxes,fontsize=10)
		ax[j].tick_params(axis="y", labelsize=10,direction="in",which='both')
		ax[j].tick_params(axis="x", labelsize=10,direction="in",which='both')
		ax[j].set_yticks([1,10,100])
		ax[j].set_ylim(1e-1,1e3)
	plt.subplots_adjust(left = 0.15,bottom = 0.1,right=0.7,top=0.9)
	ax[-1].set_xlabel('t [yr]',fontsize=10)
	ax[int(len(dirnames)/2)].set_ylabel('Q        ',fontsize=10,rotation=0)
	ax[-1].set_xlim(10,410)
	plt.subplots_adjust(left = 0.15,bottom = 0.1,right=0.7,top=0.9)
	
		
''' and so ends the useful part of this script'''		




def radial_away(a,vector,sinki,bound):
	'''line through disk l=(xsink,ysink,zsink) + scaler*(vector)
	for every cell we get the closest distance to line 'l'
	using both R_hat=(0toCell - 0toLine)      AND     (R dot vector) = 0'''
	
	origin = np.array([a.sinkx[sinki],a.sinky[sinki],a.sinkz[sinki]])
	x=a.x[bound] * code_units.d_cu
	y=a.y[bound] * code_units.d_cu
	z=a.z[bound] * code_units.d_cu
	xsink=a.sinkx[sinki] * code_units.d_cu
	ysink=a.sinky[sinki] * code_units.d_cu
	zsink=a.sinkz[sinki] * code_units.d_cu
	scaler = (vector[0]*x - vector[0]*xsink + vector[1]*y - vector[1]*ysink + vector[2]*z - vector[2]*zsink) / (vector[0]+vector[1]+vector[2])
	dx=x-xsink-scaler*vector[0]
	dy=y-ysink-scaler*vector[1]
	dz=z-zsink-scaler*vector[2]
	rs=np.sqrt(dx**2+dy**2+dz**2)

	return rs 

def define_disk(a,rho_sink,zoom):
	mask=sink_mask.sinkmask(a,zoom)[0]
	bound_to=np.ones(len(a.x))*-1
	M=a.mass*code_units.M_cu
	Msink=a.sinkmass*code_units.M_cu
	for i in range(len(a.x[mask])):
		r=np.sqrt((a.x[mask[i]]-a.sinkx)**2+(a.y[mask[i]]-a.sinky)**2+(a.z[mask[i]]-a.sinkz)**2) *code_units.d_cu
		forces=ap.G.cgs.value * M[mask[i]] * Msink  /r
		bound_to[mask[i]]=np.where(forces==forces.max())[0]


	#create space for disk data 
	sinks=[] #create enough lists to hold disc cells for all sinks
	for i in range(a.npart[-1]):
		sinks.append([])	

	rho_old=rho_sink
	for i in range(int(a.npart[-1])):
		vector=fifty_close(a,i)	#vector perpendicular to disk
		bound=np.where(bound_to==i)[0]
		rs=radial_away(a,vector,i,bound) #rs in the plane of the disk
		args=np.argsort(rs)
		keep=np.array([]) #store the disk elements here
		for j in range(len(args)):
			rho=a.rho[bound[args[j]]]
			if rho>rho_old/5:
				sinks[i].append(bound[args[j]])
				rho_old = rho




	
		







def find_disc(a,zoomzone,sink_rho):
	mask=sink_mask.sinkmask(a,zoomzone) #zooms in around the largest sink to cut down on #cells 
	mask1=np.where(a.rho[mask]>sink_rho/10) #selects highest density gas within zoomzone

	sinks=[] #create enough lists to hold disc cells for all sinks 
	for i in range(a.npart[-1]):
		sinks.append([])

	for i in range(len(a.x[mask][mask1])): #loop around the high density cells 
		f=0 #grav potential energy 
		sinkchoice='none' #belongs to no sink so far 

		Msink=a.sinkmass * code_units.M_cu
		M=a.mass[mask][mask1][i] * code_units.M_cu

		#don't allow discs to be larger than the distance to the next closest sink
		R=np.sqrt((a.x[mask][mask1][i]-a.sinkx)**2+(a.y[mask][mask1][i]-a.sinky)**2+(a.z[mask][mask1][i]-a.sinkz)**2) *code_units.d_cu

		fnew=ap.G.cgs.value * M * Msink / R
		j=np.where(fnew==fnew.max())[0][0]

		if R[j]<=np.sort(np.sqrt((a.sinkx-a.sinkx[j])**2+(a.sinky-a.sinky[j])**2+(a.sinkz-a.sinkz[j])**2))[1] *code_units.d_cu:
			sinkchoice=j

		#at this point we know which sink the cell belongs to 
		if sinkchoice != 'none':
			sinks[sinkchoice].append(mask[0][mask1[0]][i]) #put it in the right list

	plt.figure()
	arepo_utils.arepoimage(a.x[mask],a.z[mask],a.rho[mask])
	for i in range(a.npart[-1]):
		plt.scatter(a.x[sinks[i]],a.z[sinks[i]],s=0.1)
	plt.figure()
	arepo_utils.arepoimage(a.x[mask],a.y[mask],a.rho[mask])
	for i in range(a.npart[-1]):
		plt.scatter(a.x[sinks[i]],a.y[sinks[i]],s=0.1)
		




def attempt2(a,zoomzone,sink_rho):
	sinks=[]
	for i in range(a.npart[-1]):
		sinks.append([])

	mask=sink_mask.sinkmask(a,zoomzone)
	M=a.mass[mask] * code_units.M_cu
	Msink=a.sinkmass *code_units.M_cu
	plt.figure()
	for i in range(a.npart[-1]):
		print('sink '+str(i))
		r=np.sqrt((a.x[mask]-a.sinkx[i])**2+(a.y[mask]-a.sinky[i])**2+(a.z[mask]-a.sinkz[i])**2) *code_units.d_cu
		args=np.argsort(r)
		j=0
		rho_old = sink_rho
		marker=0
		u=np.array([])
		for j in range(len(a.x[mask])):
			if a.rho[mask][args[j]]>=0.6*rho_old:
				r_all=np.sqrt((a.x[mask][args[j]]-a.sinkx)**2+(a.y[mask][args[j]]-a.sinky)**2+(a.z[mask][args[j]]-a.sinkz)**2) *code_units.d_cu
				f= ap.G.cgs.value * M[args[j]] * Msink	/r_all
				if f.max()==f[i]:
					sinks[i].append(mask[0][args[j]])
					rho_old = a.rho[mask][args[j]]
				
	plt.figure()
	arepo_utils.arepoimage(a.x[mask],a.z[mask],a.rho[mask])
	for i in range(a.npart[-1]):
		plt.scatter(a.x[sinks[i]],a.z[sinks[i]],s=0.1)
	plt.figure()
	arepo_utils.arepoimage(a.x[mask],a.y[mask],a.rho[mask])
	for i in range(a.npart[-1]):
		plt.scatter(a.x[sinks[i]],a.y[sinks[i]],s=0.1)
