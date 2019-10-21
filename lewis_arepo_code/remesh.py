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


filename='/scratch/c.c1521474/rotation_collapse/ics/input_remesh.dat'                #read remesh data
a=gadget_reader_lewis.reader(filename)        

                                            #header data
sofar=[]
npart=(4000,0,0,0,0,0)
massarr=(0,0,0,0,0,0)
time=10
redshift=0
flag_sfr=0
flag_feedback=0
npartTotal=(4000,0,0,0,0,0)
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



m=a.m
v=(a.vx,a.vy,a.vz)
v_r=v_r=velocities.vary_rotation(6,(a.x,a.y,a.z),0.5,m) #m given from result of remesh
v=(v[0]+v_r[0], v[1]+v_r[1], v[2]+v_r[2])

x=(a.x,a.y,a.z)
ids=a.ids
u=a.u


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
arepo_input_writer.writer(sofar,'/scratch/c.c1521474/rotation_collapse/vary_collapse/arepo_input.dat')
