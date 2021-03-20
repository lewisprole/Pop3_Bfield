import numpy as np 
import struct

def writer(filename,bx,by,bz):
	print('writing')
	nx,ny,nz=struct.pack('i',len(bx[:,0,0])),struct.pack('i',len(bx[0,:,0])),struct.pack('i',len(bx[0,0,:]))
	sofar=nx+ny+nz
	
	blockx = struct.pack('d'*len(bx[:,0,0])**3,*bx.reshape(len(bx[:,0,0])**3))
	blocky = struct.pack('d'*len(by[0,:,0])**3,*by.reshape(len(by[0,:,0])**3))
	blockz = struct.pack('d'*len(bz[0,0,:])**3,*bz.reshape(len(bz[0,0,:])**3))
	sofar=sofar+blockx+blocky+blockz
	f= open(filename,"w+b")
	f.write(sofar)

def reader(filename):
	print('reading')
	with open(filename, mode='rb') as file:
		data = file.read()
	nx=struct.unpack('i',data[0:4])
	ny=struct.unpack('i',data[4:8])
	nz=struct.unpack('i',data[8:12])
	nx,ny,nz=nx[0],ny[0],nz[0]
	print('shape: '+str(nx)+' '+str(ny)+' '+str(nz))

	start=12
	bx=struct.unpack('d'*nx**3,data[start:start+nx**3*8])
	start=start+nx**3*8
	by=struct.unpack('d'*ny**3,data[start:start+ny**3*8])
	start=start+ny**3*8
	bz=struct.unpack('d'*nz**3,data[start:start+nz**3*8])

	bx=np.asarray(bx).reshape(nx,nx,nx)
	by=np.asarray(by).reshape(ny,ny,ny)
	bz=np.asarray(bz).reshape(nz,nz,nz)
	return bx,by,bz
