#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 12:06:08 2019

@author: lewisprole
"""

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

T=10
M=3*1.989e33
r,rho_outter=calculate_radius.BE_radius('mass',M,T)				#Bonnor Ebert radius estimation from provided mass 
boxsize=4*r							
Bsize_CU=round(boxsize/code_units.d_cu,3)
boxsize=Bsize_CU*code_units.d_cu
print('boxsize: '+ str(Bsize_CU))



x,y,z,vol_cell,vol_cell_bg=spherical_spray.uniform_sphere(1e5,1e5,r,boxsize)    #Uniform grid with +resolution sphere
print('CODE UNITS')
print('number of cells: ' +str(len(x)))
print('radius of sphere: '+str(r/code_units.d_cu))
print('cell volume insphere: '+str(vol_cell/code_units.d_cu**3))
print('density at R: '+str(rho_outter/code_units.rho_cu))
print('lowest mass in sphere: '+str(vol_cell*rho_outter/code_units.M_cu))
#x,y,z=spherical_spray.spherical_cloud(1e5,1e5,r,boxsize,boxsize,boxsize)       #Random xyz particles +resolution sphere
rho=radial_density.BE_profile(x,y,z,boxsize,T,rho_outter)			#Bonnor Ebert density profile
mid=boxsize/2
rs=np.sqrt((mid-x)**2+(mid-y)**2+(mid-z)**2)
plt.figure(),plt.plot(rs,rho,'x')
#m,rs,rho=mass.bonnor_ebert(boxsize,(x,y,z),vol_cell,vol_cell_bg,T,r)		#Or Bonnor Ebert mass profile (instead of density profile)
ids =np.linspace(1,len(x),len(x)).astype(int)
U=internal_energy.int_en(len(x),T)
v=velocities.zero_vel(len(x)) 							#0 velocities
#v=velocities.vary_rotation(boxsize,(x,y,z),0.05,m)				#Or rotation about z axis (each cell balances current G-energy within its r)  


v1=v[0]/code_units.v_cu								#convert to code units 
v2=v[1]/code_units.v_cu
v3=v[2]/code_units.v_cu
v=v1,v2,v3
#m=m/code_units.M_cu
x=x/code_units.d_cu
y=y/code_units.d_cu
z=z/code_units.d_cu
rho=rho/code_units.rho_cu





sofar=[]
npart=(len(x),0,0,0,0,0)
massarr=(0,0,0,0,0,0)
time=10
redshift=0
flag_sfr=0
flag_feedback=0
npartTotal=(len(x),0,0,0,0,0)
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

#write ICs file
sofar=arepo_input_writer.header(sofar,npart,massarr,time,redshift,flag_sfr,flag_feedback,
           npartTotal,flag_cooling,num_files,boxsize,cos1,cos2,
           hubble_param,flag_stellarage,flag_metals,npartHighword,
           flag_entropy,flag_dp,flag_1pt,scalefactor)

sofar=arepo_input_writer.tag_block(sofar,(x,y,z),'POS ','d',3)
sofar=arepo_input_writer.tag_block(sofar,v,'VEL ','d',3)
sofar=arepo_input_writer.tag_block(sofar,ids,'ID  ','i',1)
sofar=arepo_input_writer.tag_block(sofar,rho,'MASS','d',1)
sofar=arepo_input_writer.tag_block(sofar,U,'U   ','d',1)
arepo_input_writer.writer(sofar,'/scratch/c.c1521474/realistic_sphere/ics/pre_remesh.dat')
