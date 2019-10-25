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


T=10
r=int(calculate_radius.BE_radius(10*1.989e33,T))
boxsize=4*r#ap.pc.cgs.value
print('boxsize: '+ str(boxsize))
x,y,z,vol_cell,vol_cell_bg=spherical_spray.uniform_sphere(1e4,1e4,r,boxsize)
#x,y,z=spherical_spray.spherical_cloud(10000,10000,1,6,6,6)
ids =np.linspace(1,len(x),len(x)).astype(int)
U=internal_energy.int_en(len(x),T)
m,rs,rho=mass.bonnor_ebert(boxsize,(x,y,z),vol_cell,vol_cell_bg,T,r)
#v=velocities.zero_vel(len(x)) #0 velocities
v=velocities.vary_rotation(boxsize,(x,y,z),0.5,m)

v1=v[0]/code_units.v_cu
v2=v[1]/code_units.v_cu
v3=v[2]/code_units.v_cu
v=v1,v2,v3
m=m/code_units.M_cu
x=x/code_units.d_cu
y=y/code_units.d_cu
z=z/code_units.d_cu


#rho,rs=radial_density.rhos(x,y,z,6,6,6,2,2,1,-2,0)
#Mtot=radial_density.tmass(2,2,1,-2,0)


#


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
arepo_input_writer.writer(sofar,'/scratch/c.c1521474/realistic_sphere/bonnor_ebert2/arepo_input.dat')
