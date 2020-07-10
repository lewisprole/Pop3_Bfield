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
from scipy.interpolate import RegularGridInterpolator 

def turb_read(name):
	'''reading binary file
	file is in Fortran format 
	i.e |data length in bytes (int32)| -> |data (doubles)| -> |data length (int32)|
	read in this structure for x,y,z velocities'''

	with open(name, mode='rb') as file:
		data = file.read()
	buff=data[0:4] #data length in bytes 
	length=struct.unpack('i',buff)[0] #convert from binary 
	length_real=int(length/8) #data length (number of values)
	block=data[4:4+length] #read in data
	check=data[4+length:8+length] #final data length read 
	xturb=struct.unpack(length_real*'d',block) #convert data from binary 
	if buff==check: #check if the data length at the start and end are the same 
		print('x turb read')

	buff1=data[8+length:8+length+4] #repeat for y vels 
	length1=struct.unpack('i',buff1)[0]
	length1_real=int(length1/8)
	block1=data[8+length+4:8+length+4+length1]
	check1=data[8+length+4+length1:8+length+4+length1+4]
	yturb=struct.unpack(length1_real*'d',block1)
	if buff1==check1:
		print('y turb read')
	
	buff2=data[8+length+4+length1+4:8+length+4+length1+4+4] #repeat for z vels 
	length2=struct.unpack('i',buff2)[0]
	length2_real=int(length1/8)
	block2=data[8+length+4+length1+4+4:8+length+4+length1+4+4+length2]
	check2=data[8+length+4+length1+4+4+length2:8+length+4+length1+4+4+length2+4]
	zturb=struct.unpack(length2_real*'d',block2)
	if buff2==check2:
		print('z turb read')
	
	l=round(length_real**(1/3))	
	xturb=np.array(xturb).reshape(l,l,l) #convert array into data cube 
	yturb=np.array(yturb).reshape(l,l,l)	
	zturb=np.array(zturb).reshape(l,l,l)
	
	return xturb,yturb,zturb






def interpolate_turb(xturb,yturb,zturb,x,y,z,boxsize):
	'''take grid point velocities and interpolate them so that velocities can be 
	given anywhere in the cube, apply to the positions of the AREPO cells'''

	turb_size=len(xturb[:,0,0])
	side=np.linspace(0,boxsize,turb_size)
	
	xyz=np.zeros((len(x),3))
	xyz[:,0]=x
	xyz[:,1]=y
	xyz[:,2]=z
	
	fx=RegularGridInterpolator((side,side,side),xturb)
	vx=fx((x,y,z))

	fy=RegularGridInterpolator((side,side,side),yturb)
	vy=fy((x,y,z))
	
	fz=RegularGridInterpolator((side,side,side),zturb)
	vz=fz((x,y,z))

	return vx,vy,vz





def turbulence(turbulence_name,x,y,z,boxsize):
	'''read turbulent box and interpolate to cell positions'''
	xturb,yturb,zturb=turb_read(turbulence_name)
	vx,vy,vz=interpolate_turb(xturb,yturb,zturb,x,y,z,boxsize)
	return vx,vy,vz




def rescale(vx,vy,vz,alpha,M,r):
	'''uses root mean square of velocity to normalise and rescale field 
	to give desired ratio alpha=turbE/gravE'''
	v=np.sqrt(vx**2+vy**2+vz**2)
	vmean=np.mean(v)
	required_velocity=np.sqrt(2*alpha*ap.G.cgs.value*M/r)	
	factor=required_velocity/vmean #centre on 0
	v1,v2,v3=vx*factor,vy*factor,vz*factor #mean(|v|)=required velocity
	return v1,v2,v3

def rescale_from_Vrms(v_rms,vx,vy,vz):
	'''rescales field if the rms velocity is already known'''
	v=np.sqrt(vx**2+vy**2+vz**2)
	vmean=np.mean(v)
	factor=v_rms/vmean
	v1,v2,v3=vx*factor,vy*factor,vz*factor
	return v1,v2,v2

