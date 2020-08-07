import nump as np
import h5py 

def read(filename)
	hf = h5py.File(filename, 'r')
	keys=list(hf.keys())
	particle_keys=list(hf[keys[-2]])
	sink_keys=list(hf[keys[-1]])

	class a():
		pass
	
	for key in particle_keys:
		data=np.array(hf[keys[-2]][key])
		if key=='Coordinates':
			a.x,a.y,a.z=data[:,0],data[:,1],data[:,2]
		if key=='Velocities':
			a.vx,a.vy,a.vz=data[:,0],data[;,1],data[:,2]

			
			




