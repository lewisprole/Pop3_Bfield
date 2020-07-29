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


filename='/cosma7/data/dp155/dc-prol1/mid_scaling/remesh_weak/snapshot_170'                #read remesh data
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


b=0.005
bx=np.zeros_like(a.x)
bz=b*np.ones_like(a.x)
B=(bx,bx,bz)


#write ICs file 
sofar=arepo_input_writer.header(sofar,npart,massarr,time,redshift,flag_sfr,flag_feedback,
           		npartTotal,flag_cooling,num_files,boxsize,cos1,cos2,
           		hubble_param,flag_stellarage,flag_metals,npartHighword,
           		flag_entropy,flag_dp,flag_1pt,scalefactor)    
sofar=arepo_input_writer.tag_block(sofar,(a.x,a.y,a.z),'POS ','d',3)
sofar=arepo_input_writer.tag_block(sofar,(a.vx,a.vy,a.vz),'VEL ','d',3)
sofar=arepo_input_writer.tag_block(sofar,a.ids,'ID  ','i',1)
sofar=arepo_input_writer.tag_block(sofar,a.mass,'MASS','d',1)
sofar=arepo_input_writer.tag_block(sofar,a.u,'U   ','d',1)
sofar=arepo_input_writer.tag_block(sofar,B,'BFLD','d',3)
arepo_input_writer.writer(sofar,'/cosma7/data/dp155/dc-prol1/mid_scaling/weak/4node/arepo_input.dat')
#Y1e%.0f_B1e%.0f/arepo_input_Y1e%.0f_B1e%.0f.dat'%(np.log10(Gamma[j]),np.log10(Beta[j]),np.log10(Gamma[j]),np.log10(Beta[j])))

