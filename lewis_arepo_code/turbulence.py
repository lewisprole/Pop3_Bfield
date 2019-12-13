import numpy as np
import arepo_input_writer
import velocities
import spherical_spray
import radial_density
import internal_energy
import mass
import code_units
import calculate_radius
import astropy.constants as ap
import matplotlib.pyplot as plt
import plot3D
import struct

def turb_read(name):
	#a=np.fromfile(name,dtype=np.int32)

	buff = np.fromfile(name,dtype=np.int32, count=1)
	xturb = np.fromfile(name,dtype=np.double, count=64**3)
	endbuff = np.fromfile(name,dtype=np.int32, count=1)
	if buff==endbuff:
		print('x turb read')

	buff =np.fromfile(name,dtype=np.int32,count=1)
	yturb=np.fromfile(name,dtype=np.double,count=64**3)
	endbuff=np.fromfile(name,dtype=np.int32,count=1)
	if buff==endbuff:
		print('y turb read')
	
	buff=np.fromfile(name,dtype=np.int32,count=1)
	zturb=np.fromfile(name,dtype=np.double,count=64**3)
	endbuff=np.fromfile(name,dtype=np.int32,count=1)
	if buff==endbuff:
		print('z turb read')
	
	return xturb,yturb,zturb	
'''

	buffx=int(a[0]/4)
	xturb=a[1:buffx+1]
	binary=struct.pack('i'*buffx,*xturb)
	xturb=struct.unpack('d'*buffx/2,*binary)
	buffx1=int(a[buffx+1]/4)
	if buffx==buffx1:
		print('x turb read')
		print(buffx,buffx1)

	buffy=int(a[buffx+2]/4)
	yturb=a[buffx+3:buffx+2+buffy+1]
	buffy1=int(a[buffx+2+buffy+1]/4)
	if buffy==buffy1:
		print('y turb read')
		print(buffy,buffy1)	

	buffz=int(a[buffx+2+buffy+2]/4)
	zturb=a[buffx+2+buffy+3:buffx+2+buffy+2+buffz+1]
	buffz1=int(a[buffx+2+buffy+2+buffz+1]/4)
	if buffz==buffz1:
		print('z turb read')
		print(buffz,buffz1)

	shape=buffx**(1/3)
	#vx=xturb.reshape(shape,shape,shape)
	#vy=yturb.reshape(shape,shape,shape)
	#vz=zturb.reshape(shape,shape,shape)

	return xturb,yturb,zturb
'''
