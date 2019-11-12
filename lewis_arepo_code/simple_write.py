import numpy as np
import arepo_input_writer
import velocities 
import spherical_spray 
import radial_density
import internal_energy
import code_units
import astropy.constants as ap

#box parameters
T=11
M=1.989e33
r=0.016*ap.pc.cgs.value
rho_sphere=5e-18
rho_bg=rho_sphere/100
boxsize=4*r
mid=boxsize/2
G=ap.G.cgs.value

Bsize_CU=round(boxsize/code_units.d_cu,3)
print('boxsize: '+ str(Bsize_CU))
print('radius: '+str(r/code_units.d_cu))
boxsize=Bsize_CU*code_units.d_cu

#positions
N_sphere=int(1e3)
N_bg=int(1e3)
x,y,z=spherical_spray.spherical_cloud(N_sphere,N_bg,r,boxsize,boxsize,boxsize)
rs=np.sqrt((mid-x)**2+(mid-y)**2+(mid-z)**2)
x=x/code_units.d_cu
y=y/code_units.d_cu
z=z/code_units.d_cu
x=x,y,z

#mass
m1=np.ones(N_sphere) *  rho_sphere* (4/3 * np.pi * r**3 / N_sphere)     / code_units.M_cu
m2=np.ones(N_bg) *  rho_bg* (boxsize**3 - 4/3 * np.pi * r**3) / N_bg  / code_units.M_cu
m=np.append(m1,m2)


#thermal energy
U=internal_energy.int_en(len(x[0]),T)

#zero velocities pre-remesh
v=velocities.zero_vel(len(x[0]))


#ids
ids =np.linspace(1,len(x[0]),len(x[0])).astype(int)

bx,by,bz=np.zeros_like(x[0]),np.zeros_like(x[0]),np.zeros_like(x[0])
B=bx,by,bz


sofar=[]
npart=(len(x[0]),0,0,0,0,0)
massarr=(0,0,0,0,0,0)
time=10
redshift=0
flag_sfr=0
flag_feedback=0
npartTotal=(len(x[0]),0,0,0,0,0)
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
sofar=arepo_input_writer.tag_block(sofar,m,'MASS','d',1)
sofar=arepo_input_writer.tag_block(sofar,U,'U   ','d',1)
sofar=arepo_input_writer.tag_block(sofar,B,'BFLD','d',3)
arepo_input_writer.writer(sofar,'/scratch/c.c1521474/Hannebelle2/test/arepo_input.dat')

