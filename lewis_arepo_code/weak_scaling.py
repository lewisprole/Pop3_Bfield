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
import Bfield_generate
import turbulence



#box parameters - Bonnor Ebert sphere in diffuse ISM 
T=200
mu=1
#M=1.989e33
r=1.87*ap.pc.cgs.value
n0=3.7e-20/1.83
enhance=2
#n_bg=n0*enhance/100
boxsize=4*r
mid=boxsize/2
G=ap.G.cgs.value

Bsize_CU=round(boxsize/code_units.d_cu,3)
print('boxsize: '+ str(Bsize_CU))
print('radius: '+str(r/code_units.d_cu))
boxsize=Bsize_CU*code_units.d_cu

#positions
N_sphere=int(2e6)
N_bg=int(2e6)
x,y,z=spherical_spray.spherical_cloud(N_sphere,N_bg,r,boxsize,boxsize,boxsize,'no')
rs=np.sqrt(((mid-x)**2+(mid-y)**2+(mid-z)**2).astype(float))

#masses (rough)
rho_bg=2.51e-22
rho_sphere=rho_bg*100
Nsphere=len(x)/2 
Vsphere=4/3 * np.pi * r**3
Vcell=Vsphere/Nsphere
Msphere=Vcell*rho_sphere
Vbg=boxsize**3-Vsphere
Vcell_bg=Vbg/Nsphere
Mbg=Vcell_bg*rho_bg
mask=np.where(rs<r)
m=np.zeros_like(x)
m[mask]=Msphere
mask=np.where(rs>=r)
m[mask]=Mbg

#others
ids =np.linspace(1,len(x),len(x)).astype(int)
U=internal_energy.int_en(len(x),T,mu)

#magnetic field
theta=1/3 #Kolmogorov turbulence
n_bg=rho_bg
c_s=np.sqrt(ap.k_B.cgs.value*T/(ap.m_p.cgs.value))
v_rms=np.array([0.1]) * c_s
kstar,Rm_crit,box,N,n_bg,v_turb=2*np.pi/(boxsize)*101,107,boxsize,400,n_bg,v_rms
b=Bfield_generate.uniform_from_dynamo(theta,kstar,Rm_crit,box,N,n_bg,v_rms)
Bx=np.zeros_like(x) #only in z direction
By=np.zeros_like(x)
strength=b
print(strength)
Bz=np.ones_like(x)*strength/code_units.B_cu
B=(Bx,By,Bz)

#turbulence
tname='/cosma/home/dp155/dc-prol1/turbulence/vel3D.bin'
v1,v2,v3=turbulence.turbulence(tname,x,y,z,boxsize)
v1,v2,v3=turbulence.rescale_from_Vrms(v_rms,v1,v2,v3)
v1=v1/code_units.v_cu
v2=v2/code_units.v_cu
v3=v3/code_units.v_cu
v=(v1,v2,v3)

#convert to code units
x=x/code_units.d_cu
y=y/code_units.d_cu
z=z/code_units.d_cu
m=m/code_units.M_cu


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
sofar=arepo_input_writer.tag_block(sofar,m,'MASS','d',1)
sofar=arepo_input_writer.tag_block(sofar,U,'U   ','d',1)
sofar=arepo_input_writer.tag_block(sofar,B,'BFLD','d',3)
arepo_input_writer.writer(sofar,'/cosma7/data/dp155/dc-prol1/scaling_test/weak/arepo_input_'+str(len(x))+'.dat')#'/scratch/c.c1521474/popIII/Prole/ics/pre_remesh.dat')






