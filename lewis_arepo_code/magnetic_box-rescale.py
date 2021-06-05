import numpy as np
import arepo_utils 
import arepo_input_writer
import field_maker
import code_units
import read_custom_field
import matplotlib.pyplot as plt 

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

def shifter(a,box_old,zoomzone):
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
	
filename='/scratch/c.c1521474/resolution_test/MHD/1e10MHD_old/arepo_input.dat'
for W in range(1):
	a=arepo_utils.aread(filename)

	n0,n1,n2,n3,n4,n5=a.npart
	n0=len(a.x)
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
	

	for V in range(1):
		#rescale the b field strength 
		zoomzone=0.5
		Etot=sum(0.5*a.mass*(a.vx**2+a.vy**2+a.vz**2))
		B=np.sqrt(2*Etot/(2*zoomzone)**3) 
		b1,b2,b3=field_maker.rescale(a.bfield[:,0],a.bfield[:,1],a.bfield[:,2],B) #twisted field 
	
		sofar=arepo_input_writer.header(sofar,npart,massarr,time,redshift,flag_sfr,flag_feedback,
		npartTotal,flag_cooling,num_files,boxsize,cos1,cos2,
		hubble_param,flag_stellarage,flag_metals,npartHighword,
		flag_entropy,flag_dp,flag_1pt,scalefactor)
		sofar=arepo_input_writer.tag_block(sofar,(a.x,a.y,a.z),'POS ','d',3)
		sofar=arepo_input_writer.tag_block(sofar,(a.vx,a.vy,a.vz),'VEL ','d',3)
		sofar=arepo_input_writer.tag_block(sofar,a.partid,'ID  ','i',1)
		sofar=arepo_input_writer.tag_block(sofar,a.mass,'MASS','d',1)
		sofar=arepo_input_writer.tag_block(sofar,a.u,'U   ','d',1)
		sofar=arepo_input_writer.tag_block(sofar,a.rho,'RHO ','d',1)
		sofar=arepo_input_writer.tag_block(sofar,a.potential,'POT ','d',1)
		sofar=arepo_input_writer.tag_block(sofar,(a.accel[:,0],a.accel[:,1],a.accel[:,2]),'ACCE','d',3)
		sofar=arepo_input_writer.tag_block(sofar,(b1,b2,b3),'BFLD','d',3)
		sofar=arepo_input_writer.tag_block(sofar,a.divb,'DIVB','d',1)
		sofar=arepo_input_writer.tag_block(sofar,a.divbalt,'DVBA','d',1)
		sofar=arepo_input_writer.tag_block(sofar,a.psi,'BPSI','d',1)
#		sofar=arepo_input_writer.tag_block(sofar,a.dednerv,'VDED','d',1)
		#chemical species slightly trickier to prepare
		chem=np.zeros((a.chem.shape[1],len(a.chem[:,0])))
		for i in range(a.chem.shape[1]):
			chem[i]=a.chem[:,i]
		sofar=arepo_input_writer.tag_block(sofar,chem,'CHEM','d',int(a.chem.shape[1]))
		sofar=arepo_input_writer.tag_block(sofar,a.gamma,'GAMM','d',1)
		sofar=arepo_input_writer.tag_block(sofar,a.divv,'DIVV','d',1)
		sofar=arepo_input_writer.tag_block(sofar,a.peak,'PEAK','i',1)
		saveloc='/scratch/c.c1521474/resolution_test/MHD/arepo_input_weak.dat'
		arepo_input_writer.writer(sofar,saveloc)

