import numpy as np
import matplotlib.pyplot as plt
import struct
import binascii
import arepo_input_writer
import gadget_reader_lewis
import internal_energy

S=100
N=400

x=np.random.uniform(0.001,S-0.001,N)
y=np.random.uniform(0.001,S-0.001,N)
z=np.random.uniform(0.001,S-0.001,N)

rho=1e5
m=np.ones_like(x) * (rho*S**3/N)

vx,vy,vz=np.zeros_like(x),np.zeros_like(x),np.zeros_like(x)

T=800
U=internal_energy.int_en(N-2,T,1) 

ids=np.linspace(1,N,N).astype(int)

b=np.zeros(N-2)
B=(b,b,b)


#sinks
x[-1],y[-1],z[-1]=50,50,45
x[-2],y[-2],z[-2]=50,50,55

vx[-1],vy[-1],vz[-1]=0,0,1
vx[-2],vy[-2],vz[-2]=0,0,-1

m[-1],m[-2]=m[-1]*10,m[-2]*8

v=(vx,vx,vx)
X=(x,y,z)

sofar=[]
npart=(N-2,0,0,0,0,2)
massarr=(0,0,0,0,0,0)
time=10
redshift=0
flag_sfr=0
flag_feedback=0
npartTotal=(N-2,0,0,0,0,2)
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
sofar=arepo_input_writer.tag_block(sofar,X,'POS ','d',3)
sofar=arepo_input_writer.tag_block(sofar,v,'VEL ','d',3)
sofar=arepo_input_writer.tag_block(sofar,ids,'ID  ','i',1)
sofar=arepo_input_writer.tag_block(sofar,m,'MASS','d',1)
sofar=arepo_input_writer.tag_block(sofar,U,'U   ','d',1)
sofar=arepo_input_writer.tag_block(sofar,B,'BFLD','d',3)
arepo_input_writer.writer(sofar,'/scratch/c.c1521474/merger_test/arepo_input.dat')


