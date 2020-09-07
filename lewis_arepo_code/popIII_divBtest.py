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

bx,by,bz,x,y,z=field_maker.real_space(100,-2,Bsize_CU)
dx=x[1,0,0]-x[0,0,0]
x=np.linspace(dx/2,Bsize_CU-dx/2,200)
y,x,z=np.meshgrid(x,x,x)
x,y,z=x.flatten(),y.flatten(),z.flatten()
X=(x,y,z)
bx,by,bz=bx.flatten(),by.flatten(),bz.flatten()
B=(bx,by,bz)
v=velocities.zero_vel(len(x))
m=np.ones_like(x)
vx,vy,vz=v[0],v[1],v[2]
v=(vx,vy,vz)
U=internal_energy.int_en(len(x),T,1)
ids =np.linspace(1,len(x),len(x)).astype(int)


#header info
n0,n1,n2,n3,n4,n5=(len(x),0,0,0,0,0)
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

#write ics 
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
arepo_input_writer.writer(sofar,'/scratch/c.c1521474/magnetised_jeans_test/divBtest.dat')

