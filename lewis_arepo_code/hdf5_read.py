import numpy as np
import h5py 

def read(filename):
	file = h5py.File(filename,'r')
	hdf5_menu=list(file.keys())
	particle_keys=list(file['PartType0'])
	if 'PartType5' in hdf5_menu:
		sink_keys=list(file['PartType5'])

	class a():
		pass
	header=np.array(list(file.get("Header").attrs.values()))
	a.boxsize=header[6]
	npart=header[0]
	a.npart=npart
	a.nparttot=header[1]
	a.time=header[4]
	print('boxsize: '+str(a.boxsize))
	print('Npart: '+str(a.npart))
	print('Npart_tot: '+str(a.nparttot))
	print('time in code units: '+str(a.time))

	N = int(sum(npart)) - npart[4] # type 4 is reserved for TRACER_MC
	ngas = npart[0]
	nsink = npart[5]
	
	a.nsink=nsink
	a.ngas=ngas
	a.unit_leng_cm = 1.0e+17
	a.unit_mass_g = 1.991e33
	a.unit_time_s = 2.7436898e12

	igot_pot = 0
	igot_accel = 0
	igot_dt = 0
	igot_bmag = 0
	igot_chem = 0
	igot_soft = 0

	for key in particle_keys:
		data=np.array(file['PartType0'][key])
		if key=='Coordinates':
			print('Reading Coordinates')
			a.x,a.y,a.z=data[:,0],data[:,1],data[:,2]
		if key=='Velocities':
			print('Reading Velocities')
			a.vx,a.vy,a.vz=data[:,0],data[:,1],data[:,2]
		if key=='Acceleration':
			print('Reading Acceleration')
			a.accelx,a.accely,a.accelz=data[:,0],data[:,1],data[:,2]
			igot_accel = 1
		if key=='ChemicalAbundances':
			print('Reading ChemicalAbundances')
			a.chem=data			
		if key=='DednerSpeed':
			print('Reading DednerSpeed')
			a.dedner=data
		if key=='Density':
			print('Reading Density')
			a.rho=data
		if key=='Gamma':
			print('Reading Gamma')
			a.gamma=data
		if key=='InternalEnergy':
			print('Reading InternalEnergy')
			a.u=data
		if key=='MagneticField':
			print('Reading MagneticField')
			a.bx,a.by,a.bz=data[:,0],data[:,1],data[:,2]
		if key=='MagneticFieldDivergence':
			print('Reading MagneticFieldDivergence')
			a.divb=data
		if key=='MagneticFieldDivergenceAlternative':
			print('Reading MagneticFieldDivergenceAlternative')
			a.divb_alt=data
		if key=='MagneticFieldPsi':
			print('Reading MagneticFieldPsi')
			print('MagneticFieldPsi shape '+str(data.shape))
			a.bpsi=data
		if key=='Masses':
			print('Reading Masses')
			a.mass=data
		if key=='ParticleIDs':
			print('Reading ParticleIDs')
			a.partid=data
		if key=='Potential':
			print('Reading Potential')
			a.potential=data
			igot_pot = 1
		if key=='PotentialPeak':
			print('Reading PotentialPeak')
			a.peak=data
		if key=='VelocityDivergence':
			print('Reading VelocityDivergence')
			a.divv=data
	if (a.nsink > 0):
		print("Sinks read. Making sink arrays.")
		a.idsink = np.linspace(0,nsink-1,nsink,dtype='int32') + ngas
		a.sinkx = a.x[a.idsink]
		a.sinky = a.y[a.idsink]
		a.sinkz = a.z[a.idsink]
		a.sinkvx = a.vx[a.idsink]
		a.sinkvy = a.vy[a.idsink]
		a.sinkvz = a.vz[a.idsink]
		a.sinkmass = a.mass[a.idsink]
		a.sinkid = a.partid[a.idsink]

	i_not_gas = [1, 2, 3, 5] # type 4 is reservered for TRACER_MC
	if(sum(npart[i_not_gas]) > 0):
		igas = np.linspace(0,ngas-1,ngas,dtype='int32')
		a.x = a.x[igas]
		a.y = a.y[igas]
		a.z = a.z[igas]
		a.vx = a.vx[igas]
		a.vy = a.vy[igas]
		a.vz = a.vz[igas]
		a.partid = a.partid[igas]
		a.mass = a.mass[igas]
		if (igot_pot == 1):
			a.potential = a.potential[igas]
		if (igot_accel == 1):
			a.accel = a.accel[igas, :]
		if (igot_dt == 1):
			a.dt = a.dt[igas]
		if (igot_soft == 1):
			a.softening= a.softening[igas]

	# if we have chemistry, then we need to provide the temperature
	if (igot_chem > 0):
		print('Creating an array with T [K] from specific energies')
		ABHE = 0.1
		uenergy = 2.64481e+42
		ulength = 1e17
		udensity = 1.991e-18
		yn = a.rho*udensity / ((1.0 + 4.0 * ABHE) * mp)
		energy = a.u * a.rho * uenergy / ulength**3
		yntot = (1.0 + ABHE - a.chem[:, 0] + a.chem[:, 1]) * yn
		a.temp = 2.0 * energy / (3.0 * yntot * k_B)

	# get radius from densest point, com, or first sink
	if (a.nsink > 0):
		a.rad = np.sqrt((a.x - a.sinkx[0])**2 + (a.y - a.sinky[0])**2 + (a.z - a.sinkz[0])**2)


	if(igot_bmag > 0):
		a.l = (a.mass / a.rho)**(1./3.)
		a.bmag = np.sqrt(a.bx**2 + a.by[:,1]**2 + a.bz[:,2]**2)	

	return a 






