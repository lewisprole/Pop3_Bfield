import numpy as np
import h5py 

def read(filename):
	file = h5py.File('snap_000.hdf5','r')
	keys=list(file.keys())
	particle_keys=list(file[keys[3]])
	if 'PartType5' in keys:
		sink_keys=list(file[keys[4]])

	class a():
		pass
	header=np.array(list(file.get("Header").attrs.values()))
	a.boxsize=header[0]
	a.npart=header[14]
	a.time=header[20]
	for key in particle_keys:
		data=np.array(hf[keys[-2]][key])
		if key=='Coordinates':
			a.x,a.y,a.z=data[:,0],data[:,1],data[:,2]
		if key=='Velocities':
			a.vx,a.vy,a.vz=data[:,0],data[:,1],data[:,2]

			
			




