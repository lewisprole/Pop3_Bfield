import arepo_utils 
import numpy as np
import glob
import matplotlib.pyplot as plt
def sinkcheck(dirname):
	names=np.asarray(glob.glob(dirname+'/snapshot_*'))
	NAMES=[]
	for i in range(len(names)):
		N=names[i].rfind('_')
		number=names[i][N+1:]
		NAMES=np.append(NAMES,number)
	
	args=np.asarray(NAMES.argsort()).astype(int)
	print(names[args])	
	n=0
	i=0
	while n==0:
		I=args[i]
		print(names[I])
		a=arepo_utils.aread(names[I])
		n=a.nsink
		i+=1
	print('first sink in '+str(names[I]))

def rhocheck(filename):
	a=arepo_utils.aread(filename)
	if len(a.sinkx)>1:
		for i in range (len(a.sinkx)):
			core_xyz=a.sinkx[i],a.sinky[i],a.sinkz[i]
			midx=core_xyz[0]
			midy=core_xyz[1]
			midz=core_xyz[2]
			rs=np.sqrt((midx-a.x)**2+(midy-a.y)**2+(midz-a.z)**2)
			plt.figure(),plt.scatter(rs,a.rho,s=0.1)
	else:
		core_xyz=a.sinkx,a.sinky,a.sinkz
		midx=core_xyz[0]
		midy=core_xyz[1]
		midz=core_xyz[2]
		rs=np.sqrt((midx-a.x)**2+(midy-a.y)**2+(midz-a.z)**2)
		plt.figure(),plt.scatter(rs,a.rho,s=0.1)
