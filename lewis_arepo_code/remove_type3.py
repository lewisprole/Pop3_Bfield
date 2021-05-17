import numpy as np
import arepo_utils 
import arepo_input_writer

'''WARNING - you will need to keep the sink_snap file for this to work'''



filename='/scratch/c.c1521474/resolution_test/MHD2/1e10MHD/snapshot_122_original'
a=arepo_utils.aread(filename)


n0,n1,n2,n3,n4,n5=a.npart
n0=len(a.x)
n5=len(a.sinkx)
sofar=[]
npart=(n0,n1,n2,0,n4,n5)
massarr=(0,0,0,0,0,0)
time=a.time
redshift=0
flag_sfr=0
flag_feedback=0
npartTotal=(n0,n1,n2,0,n4,n5)
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

x,y,z=np.append(a.x,a.sinkx),np.append(a.y,a.sinky),np.append(a.z,a.sinkz)
vx,vy,vz=np.append(a.vx,a.sinkvx),np.append(a.vy,a.sinkvy),np.append(a.vz,a.sinkvz)
mass=np.append(a.mass,a.sinkmass)
partid=np.append(a.partid,a.sinkid)
pot=np.append(a.potential,a.sinkpotential)
accelx=np.append(a.accel[:,0],a.sinkaccelx)
accely=np.append(a.accel[:,1],a.sinkaccelx)
accelz=np.append(a.accel[:,2],a.sinkaccelz)


sofar=arepo_input_writer.header(sofar,npart,massarr,time,redshift,flag_sfr,flag_feedback,
                        npartTotal,flag_cooling,num_files,boxsize,cos1,cos2,
                        hubble_param,flag_stellarage,flag_metals,npartHighword,
                        flag_entropy,flag_dp,flag_1pt,scalefactor)
sofar=arepo_input_writer.tag_block(sofar,(x,y,z),'POS ','d',3)
sofar=arepo_input_writer.tag_block(sofar,(vx,vy,vz),'VEL ','d',3)
sofar=arepo_input_writer.tag_block(sofar,partid,'ID  ','i',1)
sofar=arepo_input_writer.tag_block(sofar,mass,'MASS','d',1)
sofar=arepo_input_writer.tag_block(sofar,a.u,'U   ','d',1)
sofar=arepo_input_writer.tag_block(sofar,a.rho,'RHO ','d',1)
sofar=arepo_input_writer.tag_block(sofar,pot,'POT ','d',1)
sofar=arepo_input_writer.tag_block(sofar,(accelx,accely,accelz),'ACCE','d',3)
sofar=arepo_input_writer.tag_block(sofar,(a.bfield[:,0],a.bfield[:,1],a.bfield[:,2]),'BFLD','d',3)
sofar=arepo_input_writer.tag_block(sofar,a.divb,'DIVB','d',1)
sofar=arepo_input_writer.tag_block(sofar,a.divbalt,'DVBA','d',1)
sofar=arepo_input_writer.tag_block(sofar,a.psi,'BPSI','d',1)
#sofar=arepo_input_writer.tag_block(sofar,a.dednerv,'VDED','d',1) #may need to add this back in - look at snapshot to know 

#chemical species slightly trickier to prepare
chem=np.zeros((a.chem.shape[1],len(a.chem[:,0])))
for i in range(a.chem.shape[1]):
	chem[i]=a.chem[:,i]
sofar=arepo_input_writer.tag_block(sofar,chem,'CHEM','d',int(a.chem.shape[1]))
sofar=arepo_input_writer.tag_block(sofar,a.gamma,'GAMM','d',1)
sofar=arepo_input_writer.tag_block(sofar,a.divv,'DIVV','d',1)
sofar=arepo_input_writer.tag_block(sofar,a.peak,'PEAK','i',1)
arepo_input_writer.writer(sofar,'/scratch/c.c1521474/resolution_test/MHD2/1e10MHD/snapshot_122')



