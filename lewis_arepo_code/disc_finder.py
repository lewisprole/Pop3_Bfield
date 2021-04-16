import numpy as np 
import matplotlib.pyplot as plt 
import sink_mask
import code_units 
import arepo_utils 

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
