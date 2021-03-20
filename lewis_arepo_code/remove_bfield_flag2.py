import numpy as np
import arepo_utils 
import arepo_input_writer
import field_maker
import code_units
import read_custom_field
import matplotlib.pyplot as plt 



filenames='/scratch/c.c1521474/resolution_test/merge/1e11_MHD/arepo_input.dat','/scratch/c.c1521474/resolution_test/merge/1e10_MHD/arepo_input.dat','/scratch/c.c1521474/resolution_test/merge/1e9_MHD/arepo_input.dat'

for W in range(len(filenames)):
	a=arepo_utils.aread(filenames[W])
	
	Bx=np.zeros_like(a.x)
	Bx,By,Bz=Bx,Bx,Bx

	n0,n1,n2,n3,n4,n5=a.npart
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


	#rescale the b field strength 
	
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
	sofar=arepo_input_writer.tag_block(sofar,(Bx,By,Bz),'BFLD','d',3)
	sofar=arepo_input_writer.tag_block(sofar,a.divb,'DIVB','d',1)
	sofar=arepo_input_writer.tag_block(sofar,a.divbalt,'DVBA','d',1)
	sofar=arepo_input_writer.tag_block(sofar,a.psi,'BPSI','d',1)
#	sofar=arepo_input_writer.tag_block(sofar,a.dednerv,'VDED','d',1)
	#chemical species slightly trickier to prepare
	chem=np.zeros((a.chem.shape[1],len(a.chem[:,0])))
	for i in range(a.chem.shape[1]):
		chem[i]=a.chem[:,i]
	sofar=arepo_input_writer.tag_block(sofar,chem,'CHEM','d',int(a.chem.shape[1]))
	sofar=arepo_input_writer.tag_block(sofar,a.gamma,'GAMM','d',1)
	sofar=arepo_input_writer.tag_block(sofar,a.divv,'DIVV','d',1)
	sofar=arepo_input_writer.tag_block(sofar,a.peak,'PEAK','i',1)
	savelocs='/scratch/c.c1521474/resolution_test/merge/1e11_noMHD/arepo_input.dat','/scratch/c.c1521474/resolution_test/merge/1e10_noMHD/arepo_input.dat','/scratch/c.c1521474/resolution_test/merge/1e9_noMHD/arepo_input.dat'
	arepo_input_writer.writer(sofar,savelocs[W])

