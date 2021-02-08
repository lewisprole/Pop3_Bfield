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
import read_custom_field
import field_maker


filename='/scratch/c.c1521474/resolution_test/remeshed.dat'               #read remesh data
a=gadget_reader_lewis.reader(filename)        
shutil.copyfile(filename,'/scratch/c.c1521474/magnetic_zooms/remeshed.dat')

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


vx,vy,vz=field_maker.create_nonscaled_field(30,-2)

#write the vel field (uniform grid) to text file incase we need it later
read_custom_field.writer('/scratch/c.c1521474/magnetic_zooms/velfield.txt',vx,vy,vz)

vx,vy,vz=field_maker.interpolate(vx,vy,vz,a.x,a.y,a.z,Bsize_CU)


for i in range (len(Gamma)):
	for j in range (len(Beta)):
		for k in range(len(Alpha)):

			#turbulence
			v1,v2,v3=turbulence.rescale(vx,vy,vz,Alpha[k],M,r)
			v1=v1/code_units.v_cu
			v2=v2/code_units.v_cu
			v3=v3/code_units.v_cu
			v=(v1,v2,v3)

			#magnetic field
			Bx=np.zeros_like(vx) #only in z direction 
			By=np.zeros_like(vy)
			Bz=np.zeros_like(vz)
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
			arepo_input_writer.writer(sofar,'/scratch/c.c1521474/magnetic_zooms/initial_collapse/arepo_input.dat')
#Y1e%.0f_B1e%.0f/arepo_input_Y1e%.0f_B1e%.0f.dat'%(np.log10(Gamma[j]),np.log10(Beta[j]),np.log10(Gamma[j]),np.log10(Beta[j])))

