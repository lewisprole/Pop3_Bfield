#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 10:22:36 2019

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
import plot3D
'''Pure molecular Hydrogen, Ideal gas, Isothermal sphere'''
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

#positions
N_sphere=int(1e4)
N_bg=int(1e4)
x,y,z=spherical_spray.spherical_cloud(N_sphere,N_bg,r,boxsize,boxsize,boxsize,'yes')
#rs=np.sqrt((mid-x)**2+(mid-y)**2+(mid-z)**2)

#densities 
rho1=rho_sphere*np.ones(N_sphere)
X,Y=mid-x[:N_sphere],mid-y[:N_sphere]
phi=np.zeros(N_sphere)
phi[np.where( abs(X)>0 )] = np.arctan( (Y/X)[np.where(abs(X)>0)] )
rho1=rho1*(1+0.5*np.cos(4*phi))
rho2=rho_bg*np.ones(N_bg)
rho=np.append(rho1,rho2)

#thermal energy 
U=U*np.ones_like(x)
#zero velocities pre-remesh, give velocities in the remesh file 
v=velocities.zero_vel(len(x)) 

#ids
ids =np.linspace(1,len(x),len(x)).astype(int)






#convert to code units 
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
arepo_input_writer.writer(sofar,'/scratch/c.c1521474/boss_bodenheimer/ics_4star/pre_remesh.dat')

