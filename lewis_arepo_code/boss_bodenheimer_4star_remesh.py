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

filename='/scratch/c.c1521474/boss_bodenheimer/ics_4star/snapshot_000'                #read remesh data
a=gadget_reader_lewis.reader(filename)        
shutil.copyfile(filename,'/scratch/c.c1521474/boss_bodenheimer/ics_4star/remeshed.dat')

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


#box parameters 
T=10
M=1.989e33
G=ap.G.cgs.value
r=3.2e16
Rgas = 8.314e7
mu=1
U=3/2 * Rgas*T/mu /code_units.t_cu**2
rho_sphere=1.44e-17
rho_bg=rho_sphere/100
boxsize=4*r
mid=boxsize/2
print("ratio = "+str(3/2 * Rgas*T/mu *M /(G*M**2/r)))

Bsize_CU=round(boxsize/code_units.d_cu,3)
print('boxsize: '+ str(Bsize_CU))
print('radius: '+str(r/code_units.d_cu))
boxsize=Bsize_CU*code_units.d_cu

#convert back to cgs
m=2*np.asarray(a.mass) #!!!!!!double the mass for double the density perturbation maxima?????
m_cgs=np.asarray(m)*code_units.M_cu  
x=np.asarray(a.x)*code_units.d_cu 
y=np.asarray(a.y)*code_units.d_cu
z=np.asarray(a.z)*code_units.d_cu

#add solid rotation
vx,vy,vz=velocities.rot_sphere(boxsize,(x,y,z),M,r,0.2)
vx=vx/code_units.v_cu
vy=vy/code_units.v_cu
vz=vz/code_units.v_cu 
v=(vx,vy,vz)

#prepare for rewrite 
x=(a.x,a.y,a.z)
ids=a.ids
u=U*np.ones_like(a.x)

B=np.zeros_like(a.x)
B=(B,B,B)


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
sofar=arepo_input_writer.tag_block(sofar,B,'BFLD','d',3)
arepo_input_writer.writer(sofar,'/scratch/c.c1521474/boss_bodenheimer/merge/4star/arepo_input.dat')

