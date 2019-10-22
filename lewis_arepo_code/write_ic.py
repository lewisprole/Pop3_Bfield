import numpy as np
import arepo_input_writer
import velocities 
import spherical_spray 
import radial_density
import internal_energy



x=spherical_spray.spherical_cloud(2000, 2000,2,6,6,6)
#
ids =np.linspace(1,4000,4000).astype(int)
#
rho,rs=radial_density.rhos(x[0],x[1],x[2],6,6,6,2,2,1,-2,0)
#Mtot=radial_density.tmass(2,2,1,-2,0)
#
U=internal_energy.int_en(4000,10)
#
v=velocities.zero_vel(4000) #0 velocities

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

#write ICs file
sofar=arepo_input_writer.header(sofar,npart,massarr,time,redshift,flag_sfr,flag_feedback,
           npartTotal,flag_cooling,num_files,boxsize,cos1,cos2,
           hubble_param,flag_stellarage,flag_metals,npartHighword,
           flag_entropy,flag_dp,flag_1pt,scalefactor)

sofar=arepo_input_writer.tag_block(sofar,x,'POS ','d',3)
sofar=arepo_input_writer.tag_block(sofar,v,'VEL ','d',3)
sofar=arepo_input_writer.tag_block(sofar,ids,'ID  ','i',1)
sofar=arepo_input_writer.tag_block(sofar,rho,'MASS','d',1)
sofar=arepo_input_writer.tag_block(sofar,U,'U   ','d',1)
arepo_input_writer.writer(sofar,'/scratch/c.c1521474/rotation_collapse/ics/pre_remesh.dat')
