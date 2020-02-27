
import numpy as np
import arepo_input_writer
import gadget_reader_lewis


Nside=50
boxsize=2
x=np.linspace(0.01,1.99,Nside)
z,y,x=np.meshgrid(x,x,x)
x=x.ravel()
y=y.ravel()
z=z.ravel()
r0=1/np.sqrt(8*np.pi)
start_mid=3*r0/2
rs=np.sqrt((start_mid-x)**2+(start_mid-y)**2)

rho=np.ones(Nside**3)
P=6
gamma=5/3
u=P/rho/(gamma-1)
cellvol=(2/Nside)**3
M=rho*cellvol
ids=np.linspace(1,Nside**3,Nside**3).astype('int')

mask=np.where(rs>r0)
Bx=1/np.sqrt(4*np.pi) * ((rs/r0)**8 - 2*(rs/r0)**4+1)
Bx[mask]=0
By=np.zeros(Nside**3)
Bz=np.ones_like(By)*1/np.sqrt(4*np.pi)

vx=np.ones_like(rho)
vy=np.ones_like(rho)
vz=np.zeros_like(rho)

#add shock
#c_s=np.sqrt(P/rho)
#c_s_alfven=Bz[0]/np.sqrt(rho[0])
#shock_mask=np.where(y<0.5)
#vy[shock_mask]=2*c_s


v=(vx,vy,vz)
x=(x,y,z)
B=(Bx,By,Bz)


#Header info 

sofar=[]
npart=(len(x[0]),0,0,0,0,0)
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

#write file

sofar=arepo_input_writer.header(sofar,npart,massarr,time,redshift,flag_sfr,flag_feedback,
               npartTotal,flag_cooling,num_files,boxsize,cos1,cos2,
               hubble_param,flag_stellarage,flag_metals,npartHighword,
               flag_entropy,flag_dp,flag_1pt,scalefactor)

sofar=arepo_input_writer.tag_block(sofar,x,'POS ','d',3)
sofar=arepo_input_writer.tag_block(sofar,v,'VEL ','d',3)
sofar=arepo_input_writer.tag_block(sofar,ids,'ID  ','i',1)
sofar=arepo_input_writer.tag_block(sofar,M,'MASS','d',1)
sofar=arepo_input_writer.tag_block(sofar,u,'U   ','d',1)
sofar=arepo_input_writer.tag_block(sofar,B,'BFLD','d',3)
arepo_input_writer.writer(sofar,'/scratch/c.c1521474/Tricco/noshock/arepo_input.dat')




