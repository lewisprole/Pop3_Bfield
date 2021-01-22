import numpy as np
import arepo_utils 
import arepo_input_writer


def masker(name,zoomzone):
	'''creates a mask giving the cells within a cube of length 2*zoomzone, to be used as 
	xnew=a.x[mask] etc'''
	a=arepo_utils.aread(name)
	print('original no cells: '+str(len(a.x)))
	mid=np.where(a.rho==a.rho.max())
	maskx=np.where(abs(a.x-a.x[mid])<zoomzone)
	masky=np.where(abs(a.y-a.y[mid])<zoomzone)
	maskz=np.where(abs(a.z-a.z[mid])<zoomzone)
	MASK=np.intersect1d(maskx,masky)
	MASK=np.intersect1d(MASK,maskz)
	print('new boxsize: '+str(2*zoomzone))
	print('no cells: '+str(len(MASK)))
	return MASK 

def shifter(a,mask,box_old,zoomzone):
	mid=np.where(a.rho==a.rho.max())
	x_to_side=a.x[mid]
	y_to_side=a.y[mid]
	z_to_side=a.z[mid]
	xshift = -x_to_side + zoomzone
	yshift = -y_to_side + zoomzone
	zshift = -z_to_side + zoomzone
	xnew=a.x+xshift 
	ynew=a.y+yshift 
	znew=a.z+zshift
	return xnew,ynew,znew 
	
	
	
#filename='/cosma7/data/dp155/dc-prol1/uniform_field/0.1/snapshot_255'
#filename='/cosma7/data/dp155/dc-prol1/resolution_test/1e8_sink/snapshot_265'
#filename='/scratch/c.c1521474/resolution_test/merge/1e8_redo/snapshot_042'
filename='/scratch/c.c1521474/resolution_test/merge/1e11/snapshot_037'
a=arepo_utils.aread(filename)
#zoomzone=0.1 
#zoomzone=1
#zoomzone=0.4
zoomzone=0.05
box=230.809
mask=masker(filename,zoomzone)
x,y,z=shifter(a,mask,box,zoomzone)
print(min(x[mask]),max(x[mask]),min(y[mask]),max(y[mask]),min(z[mask]),max(z[mask]))


n0,n1,n2,n3,n4,n5=a.npart
n0=len(mask)
sofar=[]
npart=(n0,n1,n2,n3,n4,n5)
massarr=(0,0,0,0,0,0)
time=a.time
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


sofar=arepo_input_writer.header(sofar,npart,massarr,time,redshift,flag_sfr,flag_feedback,
                        npartTotal,flag_cooling,num_files,boxsize,cos1,cos2,
                        hubble_param,flag_stellarage,flag_metals,npartHighword,
                        flag_entropy,flag_dp,flag_1pt,scalefactor)
sofar=arepo_input_writer.tag_block(sofar,(x[mask],y[mask],z[mask]),'POS ','d',3)
sofar=arepo_input_writer.tag_block(sofar,(a.vx[mask],a.vy[mask],a.vz[mask]),'VEL ','d',3)
sofar=arepo_input_writer.tag_block(sofar,a.partid[mask],'ID  ','i',1)
sofar=arepo_input_writer.tag_block(sofar,a.mass[mask],'MASS','d',1)
sofar=arepo_input_writer.tag_block(sofar,a.u[mask],'U   ','d',1)
sofar=arepo_input_writer.tag_block(sofar,a.rho[mask],'RHO ','d',1)
sofar=arepo_input_writer.tag_block(sofar,a.potential[mask],'POT ','d',1)
sofar=arepo_input_writer.tag_block(sofar,(a.accel[:,0][mask],a.accel[:,1][mask],a.accel[:,2][mask]),'ACCE','d',3)
sofar=arepo_input_writer.tag_block(sofar,(a.bfield[:,0][mask],a.bfield[:,1][mask],a.bfield[:,2][mask]),'BFLD','d',3)
sofar=arepo_input_writer.tag_block(sofar,a.divb[mask],'DIVB','d',1)
sofar=arepo_input_writer.tag_block(sofar,a.divbalt[mask],'DVBA','d',1)
sofar=arepo_input_writer.tag_block(sofar,a.psi[mask],'BPSI','d',1)
sofar=arepo_input_writer.tag_block(sofar,a.dednerv[mask],'VDED','d',1)
#chemical species slightly trickier to prepare
chem=np.zeros((a.chem.shape[1],len(a.chem[:,0][mask])))
for i in range(a.chem.shape[1]):
	chem[i]=a.chem[:,i][mask]
sofar=arepo_input_writer.tag_block(sofar,chem,'CHEM','d',int(a.chem.shape[1]))
sofar=arepo_input_writer.tag_block(sofar,a.gamma[mask],'GAMM','d',1)
sofar=arepo_input_writer.tag_block(sofar,a.divv[mask],'DIVV','d',1)

#arepo_input_writer.writer(sofar,'/cosma7/data/dp155/dc-prol1/crop_scaling/strong/4node/arepo_input.dat')
#arepo_input_writer.writer(sofar,'/cosma7/data/dp155/dc-prol1/periodic_box/arepo_input.dat')
#arepo_input_writer.writer(sofar,'/scratch/c.c1521474/resolution_test/merge/1e8_cutout/arepo_input.dat')
arepo_input_writer.writer(sofar,'/scratch/c.c1521474/resolution_test/merge/1e11_cutout/arepo_input.dat')

