import numpy as np 


def writer(filename,bx,by,bz):
	print('writing')
	array=np.hstack(( bx.ravel(),by.ravel(),bz.ravel() )).ravel()
	np.savetxt(filename, array, fmt="%s")

def reader(filename,N):
	print('reading')
	N=N*2
	data=np.loadtxt(filename)
	vx=data[:N**3].reshape(N,N,N)
	vy=data[N**3:2*N**3].reshape(N,N,N)
	vz=data[2*N**3:3*N**3].reshape(N,N,N)
	return vx,vy,vz
