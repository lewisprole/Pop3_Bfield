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
import field_maker

'''most of this doesn't matter, just trying to get an increasingly resolved random spray of cells'''

#box parameters - Bonnor Ebert sphere in diffuse ISM 
T=200
mu=1
#M=1.989e33
r=1
boxsize=4*r
mid=boxsize/2
G=ap.G.cgs.value


bx,by,bz=field_maker.create_nonscaled_Bfield(50,3/2)#10,3/2)
B=np.sqrt(bx**2+by**2+bz**2)
B=B.max()
bx,by,bz=bx/B,by/B,bz/B

vx,vy,vz=field_maker.create_nonscaled_Bfield(50,-2)#10,-2)
V=np.sqrt(vx**2+vy**2+vz**2)
V=V.max()
vx,vy,vz=vx/V,vy/V,vz/V


#positions
Nfield=50**3#10**3 
N=int(Nfield) * np.array([5e-1])#1e-1,1,10])
locs='1e-1','1','10'
for i in range(len(N)):
	print(N[i])
	x,y,z=spherical_spray.uniform_spray(int(N[i]),boxsize,boxsize,boxsize)
	Bx,By,Bz=field_maker.interpolate(bx,by,bz,x,y,z,boxsize)
	Vx,Vy,Vz=field_maker.interpolate(vx,vy,vz,x,y,z,boxsize)
	rho=np.ones_like(x)*boxsize**3/N[i]

	#others
	ids =np.linspace(1,len(x),len(x)).astype(int)
	U=internal_energy.int_en(len(x),T,mu)


	sofar=[]
	npart=(len(x),0,0,0,0,0)
	massarr=(0,0,0,0,0,0)
	time=0
	redshift=0
	flag_sfr=0
	flag_feedback=0
	npartTotal=(len(x),0,0,0,0,0)
	flag_cooling=0
	num_files=1
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
	sofar=arepo_input_writer.tag_block(sofar,(Vx,Vy,Vz),'VEL ','d',3)
	sofar=arepo_input_writer.tag_block(sofar,ids,'ID  ','i',1)
	sofar=arepo_input_writer.tag_block(sofar,rho,'MASS','d',1)
	sofar=arepo_input_writer.tag_block(sofar,U,'U   ','d',1)
	sofar=arepo_input_writer.tag_block(sofar,(Bx,By,Bz),'BFLD','d',3)
	print('/scratch/c.c1521474/interpolation_trial/'+locs[i]+'/arepo_input.dat')
	arepo_input_writer.writer(sofar,'/scratch/c.c1521474/interpolation_trial/5e-1/arepo_input.dat')






