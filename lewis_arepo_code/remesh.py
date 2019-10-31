#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 13:44:03 2019

@author: lewisprole
"""

'''script to read remesh file and extract the mass to calculate angular velocities'''

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
filename='/scratch/c.c1521474/realistic_sphere/ics/snapshot_023'                #read remesh data
a=gadget_reader_lewis.reader(filename)        

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




T=10
M=3*1.989e33
r,rho_outter=calculate_radius.BE_radius('mass',M,T)
xsize=4*r#ap.pc.cgs.value
Bsize_CU=round(xsize/code_units.d_cu,3)
xsize=Bsize_CU*code_units.d_cu
print('boxsize: '+str(Bsize_CU))



v=(a.vx,a.vy,a.vz)
#v_r=v_r=velocities.vary_rotation(r,(a.x,a.y,a.z),0.05,m) #m given from result of remesh

#convert back to cgs
m=a.mass
m_cgs=np.asarray(m)*code_units.M_cu  
x=np.asarray(a.x)*code_units.d_cu 
y=np.asarray(a.y)*code_units.d_cu
z=np.asarray(a.z)*code_units.d_cu

#add bulk velocity 
v=(a.vx,a.vy,a.vz)
vx,vy,vz=velocities.rot_sphere(xsize,(x,y,z),M,r,0.05)

#convert into code units
vx=vx/code_units.v_cu
vy=vy/code_units.v_cu
vz=vz/code_units.v_cu 
v=(v[0]+vx, v[1]+vy, v[2]+vz)

#prepare for rewrite 
x=(a.x,a.y,a.z)
ids=a.ids
u=a.u

#estimate freefall time in code units 
t_ff=calculate_freefall.ff(m,x,r/code_units.d_cu,Bsize_CU)
print('free-fall time :' +str(t_ff))

#write ICs file 
sofar=arepo_input_writer.header(sofar,npart,massarr,time,redshift,flag_sfr,flag_feedback,
           npartTotal,flag_cooling,num_files,boxsize,cos1,cos2,
           hubble_param,flag_stellarage,flag_metals,npartHighword,
           flag_entropy,flag_dp,flag_1pt,scalefactor)

sofar=arepo_input_writer.tag_block(sofar,x,'POS ','d',3)
sofar=arepo_input_writer.tag_block(sofar,v,'VEL ','d',3)
sofar=arepo_input_writer.tag_block(sofar,ids,'ID  ','i',1)
sofar=arepo_input_writer.tag_block(sofar,m,'MASS','d',1)
sofar=arepo_input_writer.tag_block(sofar,u,'U   ','d',1)
arepo_input_writer.writer(sofar,'/scratch/c.c1521474/realistic_sphere/bonnor_ebert2/arepo_input.dat')
