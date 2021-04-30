import numpy as np 
import matplotlib.pyplot as plt 
import sink_mask
import code_units 
import arepo_utils


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



def fifty_close(a,accradius,sinki):
	'''look at 50 cells closest to sink
	assume that cells in disk have higher ratio of rotational velocity to total velocity 
	take the unit vector of the rotational velocity to be a rod through the disk
	return the rotational velocity unit vector
	'''
	r=np.sqrt((a.x-a.sinkx[sinki])**2+(a.y-a.sinky[sinki])**2+(a.z-a.sinkz[sinki])**2) *code_units.d_cu
	mask=np.where(r<accradius*2 *code_units.d_cu)[0] #new
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
		rvec=np.array([ a.x[I]-a.sinkx[sinki], a.y[I]-a.sinky[sinki], a.z[I]-a.sinkz[sinki] ]) *code_units.d_cu
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

def get_disc(a,rho_sink,accradius,zoomzone):

	sinks=[] #create array for each of the sinks, holds index of cells in their disc 
	for i in range(a.npart[-1]):
		sinks.append([])

	bound_to=np.ones(len(a.x)) * -1 #array of -1's, to be set to the sinki of the sink they are bound to
	mask=sink_mask.sinkmask(a,zoomzone)[0] #only do this in the very central region of the simulation box 
	M=a.mass * code_units.M_cu
	Msink=a.sinkmass *code_units.M_cu
	for i in range(len(a.x[mask])): #loop through cells to find who they are bound to 
		r=np.sqrt((a.x[mask[i]]-a.sinkx)**2+(a.y[mask[i]]-a.sinky)**2+(a.z[mask[i]]-a.sinkz)**2) *code_units.d_cu #to all sinks
		forces=ap.G.cgs.value * M[mask[i]] * Msink  /r #grav force to all sinks
		vmag=np.sqrt(a.vx[mask[i]]**2+a.vy[mask[i]]**2+a.vz[mask[i]]**2) *code_units.v_cu
		sink_index=np.where(forces==forces.max())[0]
		if 0.5 * M[mask[i]] *vmag**2 < forces.max(): #if it's bound,check it's in circular motion
			rvec=np.array([(a.x[mask[i]]-a.sinkx[sink_index]),(a.y[mask[i]]-a.sinky[sink_index]),(a.z[mask[i]]-a.sinkz[sink_index])]) *code_units.d_cu
			vvec=np.array([a.x[mask[i]],a.y[mask[i]],a.z[mask[i]]]) *code_units.v_cu
			cross1 = vvec[1] * rvec[2] - vvec[2] * rvec[1]
			cross2 = -vvec[0] * rvec[2] + vvec[2] * rvec[0]
			cross3 = vvec[0] * rvec[1] - vvec[1] * rvec[0]
			vrot = 1/r[sink_index] * np.sqrt(cross1**2+cross2**2+cross3**2) #velocity perpendicular to r (to sink)
			if vrot>0.5*vmag: #more than half the velocity is circular	
				bound_to[mask[i]]=sink_index #index of the sink they are bound to

	for i in range(int(a.npart[-1])): #loop through each sink
		vector=fifty_close(a,accradius,i) #find vector through the disc 
		bound=np.where(bound_to==i)[0]
		dx=(a.x[bound]-a.sinkx[i]) *code_units.d_cu
		dy=(a.y[bound]-a.sinky[i]) *code_units.d_cu
		dz=(a.z[bound]-a.sinkz[i]) *code_units.d_cu
		r=np.sqrt(dx**2+dy**2+dz**2)
		vecmag=np.sqrt(vector[0]**2 + vector[1]**2 +vector[2]**2)
		costheta = (vector[0]*dx +  vector[1]*dy +  vector[2]*dz) / r / vecmag
		angle=np.where(abs(costheta)<0.5)[0]
		
		args=np.argsort(r[angle])
		rho_old=rho_sink
		marker=0
		j=0
		Qs=np.array([])
		while marker==0:
		#for j in range(len(args)):
			if j<len(args):
				if j==0:
					rho_av=a.rho[bound][angle][args[j]]
					sinks[i].append(bound[angle][args[j]])
				
				else:
					if r[angle][args[j]]<0.001*code_units.d_cu:
						if a.rho[bound][angle][args[j]]>rho_av/10:
							sinks[i].append(bound[angle][args[j]])
							#plt.plot(r[angle][args[j]],a.rho[bound][angle][args[j]],'.',c='k')
							rho_av=np.mean(a.rho[sinks[i]])
							#plt.plot(r[angle][args[j]],rho_av,'.',c='k')
							c_s=np.sqrt(ap.k_B.cgs.value * a.temp[bound][angle][args[j]]/ap.m_p.cgs.value)
							surface_density = a.mass[bound][angle][args[j]]*code_units.M_cu / ((a.mass[bound][angle][args[j]]/a.rho[bound][angle][args[j]])**(1/3)*code_units.d_cu)**2
							#units get weird, r needs to be AU, G base units, M in solar masses, gives time in years 
							period = np.sqrt(4*np.pi**2 * (r[angle][args[j]]/ap.au.cgs.value)**3 / (ap.G.value*a.sinkmass[i]*code_units.M_cu/ap.M_sun.cgs.value)) * (60*60*24*365)
							frequency = 2*np.pi/period 
							Q=c_s*frequency/(np.pi*ap.G.cgs.value*surface_density)
							Qs=np.append(Qs,Q)
							#print(len(sinks[i]),len(Qs))
						if a.rho[bound][angle][args[j]] < rho_av/10000:
							marker=1
							print('density break')
				j+=1
			else:
				marker=1
					
		R=np.sqrt((a.x[sinks[i]]-a.sinkx[i])**2+(a.y[sinks[i]]-a.sinky[i])**2+(a.z[sinks[i]]-a.sinkz[i])**2)
		plt.plot(R[1:],Qs,'.')







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
