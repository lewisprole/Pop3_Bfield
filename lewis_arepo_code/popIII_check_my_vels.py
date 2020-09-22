#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 12:07:05 2019

@author: lewisprole
"""
import numpy as np
import matplotlib.pyplot as plt 
import struct
import binascii
import arepo_input_writer
import gadget_reader_lewis
import spherical_spray
import velocities
import radial_density
import internal_energy
import calculate_radius
import code_units
import calculate_freefall
import os, shutil
import astropy.constants as ap 
import turbulence
import Bfield
import field_maker


filename='/scratch/c.c1521474/popIII/Prole/ics_turb_highres/snapshot_000'                #read remesh data
a=gadget_reader_lewis.reader(filename)        
#shutil.copyfile(filename,'/scratch/c.c1521474/popIII/Prole/resolution_test/remeshed.dat')

n0,n1,n2,n3,n4,n5=a.npart
                                            #header data
sofar=[]
npart=(n0,n1,n2,n3,n4,n5)
massarr=(0,0,0,0,0,0)
time=10
redshift=0
flag_sfr=0
flag_feedback=0
npartTotal=(n0,n1,n2,n3,n4,n5)
flag_cooling=0
num_files=1
boxsize=0
cos1=0
cos2=0
hubble_param=1
flag_stellarage=0
flag_metals=0
npartHighword=(0,0,0,0,0,0)
flag_entropy=0
flag_dp=1
flag_1pt=0
scalefactor=1


#box parameters - Bonnor Ebert sphere in diffuse ISM
T=200
r=1.87*ap.pc.cgs.value
print('Radius: '+str(r))
n0=3.7e-20/1.83
enhance=1.83
n_bg=n0*enhance/100
boxsize=4*r
mid=boxsize/2
G=ap.G.cgs.value

Bsize_CU=round(boxsize/code_units.d_cu,3)
print('boxsize: '+ str(Bsize_CU))
print('radius: '+str(r/code_units.d_cu))
boxsize=Bsize_CU*code_units.d_cu

#calculate mass in sphere 
mid=Bsize_CU/2
rs=np.sqrt((mid-a.x)**2+(mid-a.y)**2+(mid-a.z)**2)
mask=np.where(rs<r/code_units.d_cu)
M=sum(np.array(a.mass)[mask])*code_units.M_cu
print('Mass: '+str(M))

#convert back to cgs
m=a.mass
m_cgs=np.asarray(m)*code_units.M_cu  
x=np.asarray(a.x)*code_units.d_cu 
y=np.asarray(a.y)*code_units.d_cu
z=np.asarray(a.z)*code_units.d_cu
U=internal_energy.int_en(len(x),T,1)

#rotation
#Beta=np.array([1e-5,1e-4])
Beta=np.array([0]) 

#magnetic field
#Gamma=np.array([1e-2,1e-8])
Gamma=np.array([0]) 

#turbulence
Alpha=np.array([0.05])
#Alpha=np.array([0])

for i in range (len(Gamma)):
	for j in range (len(Beta)):
		for k in range(len(Alpha)):
			#add solid rotation
			if Beta[j]==0:
				v=velocities.zero_vel(len(x))
				vx,vy,vz=v[0],v[1],v[2]
			else:
				vx,vy,vz=velocities.rot_sphere(boxsize,(x,y,z),M,r,Beta[j])
				vx=vx/code_units.v_cu
				vy=vy/code_units.v_cu
				vz=vz/code_units.v_cu 
				v=(vx,vy,vz)

			#turbulence
			if Alpha[k]>0:
				#tname='/home/c.c1521474/turbulent_box/64_cube/nat_seed_43/vel3D.bin'
				#v1,v2,v3=turbulence.turbulence(tname,x,y,z,boxsize)
				#v1,v2,v3=turbulence.rescale(v1,v2,v3,Alpha[k],M,r)
				
				v1,v2,v3=field_maker.prepare_ICs(100,-2,boxsize,x,y,z,np.sqrt(2*Alpha[i]*ap.G.cgs.value*M/r))
					
			
				v1=v1/code_units.v_cu
				v2=v2/code_units.v_cu
				v3=v3/code_units.v_cu
				v=(vx+v1,vy+v2,vx+v3)
				vmag=np.sqrt( (vx+v1)**2 + (vy+v2)**2 + (vx+v3)**2)

			#magnetic field
			kenetic=sum(0.5*m_cgs*(vmag*code_units.v_cu)**2)
			magnetic=0.1*kenetic #10%
			strength=np.sqrt(2 * magnetic / boxsize**3)
			print('B field: '+str(strength)+'G')
			Bx,By,Bz=field_maker.prepare_ICs(100,3/2,boxsize,x,y,z,strength)
			Bx,By,Bz= Bx/code_units.B_cu, By/code_units.B_cu, Bz/code_units.B_cu
			B=(Bx,By,Bz)	
	

			#prepare for rewrite
			X=(a.x,a.y,a.z)
			ids=a.ids
			


#crit_MtoF=0.53/(3*np.pi) * np.sqrt(5/G) #critical mass-to-flux (cgs)
#mu=1 #ratio of mass-to-flux over critical mass-to-flux
#MtoF=mu*crit_MtoF #mass-to-flux ratio (cgs)
#F=M/MtoF #flux (cgs)
#B=F/(np.pi*(r)**2) #flux densiity (G)    
#Bcode=B/code_units.B_cu #into code units 
#Bz=Bcode*np.ones_like(vx) #only in z direction 
#B=(Bx,By,Bz)


			#write ICs file 
			sofar=arepo_input_writer.header(sofar,npart,massarr,time,redshift,flag_sfr,flag_feedback,
           		npartTotal,flag_cooling,num_files,boxsize,cos1,cos2,
           		hubble_param,flag_stellarage,flag_metals,npartHighword,
           		flag_entropy,flag_dp,flag_1pt,scalefactor)    
			sofar=arepo_input_writer.tag_block(sofar,X,'POS ','d',3)
			sofar=arepo_input_writer.tag_block(sofar,v,'VEL ','d',3)
			sofar=arepo_input_writer.tag_block(sofar,ids,'ID  ','i',1)
			sofar=arepo_input_writer.tag_block(sofar,m,'MASS','d',1)
			sofar=arepo_input_writer.tag_block(sofar,U,'U   ','d',1)
			sofar=arepo_input_writer.tag_block(sofar,B,'BFLD','d',3)
			arepo_input_writer.writer(sofar,'/scratch/c.c1521474/test_field/arepo_input.dat')

