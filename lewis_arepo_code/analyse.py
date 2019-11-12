import arepo_utils 
import numpy as np
import glob
import matplotlib.pyplot as plt
import code_units 
import astropy.constants as ap
from scipy.stats import binned_statistic
import io
import sys

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

def rhocheck(filename,boxsize):
	a=arepo_utils.aread(filename)
	rho=a.rho*code_units.rho_cu
	try:
		if len(a.sinkx)>1:
			for i in range (len(a.sinkx)):
				core_xyz=a.sinkx[i],a.sinky[i],a.sinkz[i]
				midx=core_xyz[0]
				midy=core_xyz[1]
				midz=core_xyz[2]
				rs=np.sqrt((midx-a.x)**2+(midy-a.y)**2+(midz-a.z)**2)
				rs=rs*code_units.d_cu
				rs=rs/ap.pc.cgs.value
				R=boxsize/4
				mask=np.where(np.absolute(midz-a.z)<0.1*R)
				plt.figure(),plt.scatter(np.log10(rs[mask]),np.log10(rho[mask]),s=0.1)
				plt.xlabel('log10(r)'),plt.ylabel('log10(rho)')
		else:
			core_xyz=a.sinkx,a.sinky,a.sinkz
			midx=core_xyz[0]
			midy=core_xyz[1]
			midz=core_xyz[2]
			rs=np.sqrt((midx-a.x)**2+(midy-a.y)**2+(midz-a.z)**2)
			rs=rs*code_units.d_cu
			rs=rs/ap.pc.cgs.value
			R=boxsize/4
			mask=np.where(np.absolute(midz-a.z)<0.1*R)
			plt.figure(),plt.scatter(np.log10(rs[mask]),np.log10(rho[mask]),s=0.1)
			plt.xlabel('log10(r)'),plt.ylabel('log10(rho)')
			bin_means = binned_statistic(np.log10(rs[mask]), np.log10(rho[mask]), bins=200, range=(-4.5,1.5))
			plt.plot(bin_means[1][:-1],bin_means[0],c='r')
	except AttributeError:
		
		midx=boxsize/2
		midy=boxsize/2
		midz=boxsize/2
		rs=np.sqrt((midx-a.x)**2+(midy-a.y)**2+(midz-a.z)**2)
		rs=rs*code_units.d_cu
		rs=rs/ap.pc.cgs.value
		R=boxsize/4
		mask=np.where(np.absolute(midz-a.z)<0.1*R)
		plt.figure(),plt.scatter(np.log10(rs[mask]),np.log10(rho[mask]),s=0.1)
		plt.xlabel('log10(r)'),plt.ylabel('log10(rho)')
		bin_means = binned_statistic(np.log10(rs[mask]), np.log10(rho[mask]), bins=200, range=(-4.5,1.5))
		plt.plot(bin_means[1][:-1],bin_means[0],c='r')
		#plt.plot(bin_means[1][:-1],bin_means[0],'x')

def timecheck(dirname,target):
	names=np.asarray(glob.glob(dirname+'/snapshot_*'))
	NAMES=[]
	for i in range(len(names)):
		N=names[i].rfind('_')
		number=names[i][N+1:]
		NAMES=np.append(NAMES,number)

	args=np.asarray(NAMES.argsort()).astype(int)
	output=[]
	output1=[]
	
	for I in range(len(target)):
		DIF=1e20
		text_trap = io.StringIO()
		sys.stdout = text_trap
		c=names[args[0]]
		for arg in args:
			a=arepo_utils.aread(names[arg])
			dif=np.absolute(target[I]-a.time*code_units.t_cu)
			if dif>DIF:
				output=np.append(output,c)
				output1=np.append(output1,a.time*code_units.t_cu)
				break
			DIF=dif
			c=names[arg]
		sys.stdout = sys.__stdout__
	return output,output1
		

def cylindrical_velocity(filename,boxsize):
	'''radial and rotational velocity in the xy plane at z=0 (within 0.1*R)'''
	a=arepo_utils.aread(filename)
	R=boxsize/2
	try:
		if len(a.sinkx)>1:
			for i in range (len(a.sinkx)):
				midx,midy,midz=a.sinkx[i],a.sinky[i],a.sinkz[i]
				x,y=(a.x-midx),(a.y-midy)
				rs=np.sqrt(x**2+y**2)
				vr=(x*a.vx+y*a.vy)/rs
				vr=vr*code_units.v_cu / 1e5 #code units to cgs to km/s
				vtheta=(x*a.vy-y*a.vx)/(x**2+y**2)
				vtheta=vtheta*code_units.v_cu / 1e5
				mask=np.where(np.absolute(midz-a.z)<0.1*R)
				plt.figure(),plt.scatter(np.log10(rs[mask]),vr[mask],s=0.1)#,plt.ylim(-18.5,-11.5),plt.xlim(-4.5,-1.5)
				plt.xlabel('log10(r)'),plt.ylabel('vr')
				plt.figure(),plt.scatter(np.log10(rs[mask]),vtheta[mask],s=0.1)#,plt.ylim(-18.5,-11.5),plt.xlim(-4.5,-1.5)
				plt.xlabel('log10(r)'),plt.ylabel('vtheta')
		else:
			midx,midy,midz=a.sinkx,a.sinky,a.sinkz
			x,y=(a.x-midx),(a.y-midy)
			rs=np.sqrt(x**2+y**2)
			vr=(x*a.vx+y*a.vy)/rs
			vr=vr*code_units.v_cu / 1e5 #code units to cgs to km/s
			vtheta=(x*a.vy-y*a.vx)/(x**2+y**2)
			vtheta=vtheta*code_units.v_cu / 1e5
			mask=np.where(np.absolute(midz-a.z)<0.1*R)
			plt.figure(),plt.scatter(np.log10(rs[mask]),vr[mask],s=0.1)#,plt.ylim(-18.5,-11.5),plt.xlim(-4.5,-1.5)
			plt.xlabel('log10(r)'),plt.ylabel('vr')
			plt.figure(),plt.scatter(np.log10(rs[mask]),vtheta[mask],s=0.1)#,plt.ylim(-18.5,-11.5),plt.xlim(-4.5,-1.5)
			plt.xlabel('log10(r)'),plt.ylabel('vtheta')

	except AttributeError:
		mid=boxsize/2
		x,y=(a.x-mid),(a.y-mid)
		rs=np.sqrt(x**2+y**2)
		vr=(x*a.vx+y*a.vy)/rs
		vr=vr*code_units.v_cu / 1e5 #code units to cgs to km/s
		vtheta=(x*a.vy-y*a.vx)/(x**2+y**2)
		print((vr))
		vtheta=vtheta*code_units.v_cu / 1e5
		mask=np.where(np.absolute(mid-a.z)<0.1*R)
		plt.figure(),plt.scatter(np.log10(rs[mask]),vr[mask],s=0.1)#,plt.ylim(-1.4,0),plt.xlim(-4.5,-1.5)
		plt.xlabel('log10(r)'),plt.ylabel('vr')
		plt.figure(),plt.scatter(np.log10(rs[mask]),vtheta[mask],s=0.1)#,plt.ylim(0,2),plt.xlim(-4.5,-1.5)
		plt.xlabel('log10(r)'),plt.ylabel('vtheta')
	
def Bcheck(filename,size):
	a=arepo_utils.aread(filename)
	mid=size/2
	R=size/4
	try:
		if len(a.sinkx)>1:
			midx,midy,midz=a.sinkx[i],a.sinky[i],a.sinkz[i]
			x,y=(a.x-midx),(a.y-midy)
			rs=np.sqrt(x**2+y**2)
			mask=np.where(np.absolute(midz-a.z)<0.1*R)
			B=a.bfield[:,2] * code_units.B_cu *1e6
			plt.figure()
			plt.scatter(np.log10(rs[mask]),np.log10(B[mask]),s=0.1)
			plt.xlabel('log10(r)'),plt.ylabel('log10(B_z)')
		
		else:
			midx,midy,midz=a.sinkx,a.sinky,a.sinkz
			x,y=(a.x-midx),(a.y-midy)
			rs=np.sqrt(x**2+y**2)
			mask=np.where(np.absolute(midz-a.z)<0.1*R)
			B=a.bfield[:,2] * code_units.B_cu *1e6
			plt.figure()
			plt.scatter(np.log10(rs[mask]),np.log10(B[mask]),s=0.1)
			plt.xlabel('log10(r)'),plt.ylabel('log10(B_z)')
	except AttributeError:
		x,y=(a.x-mid),(a.y-mid)
		rs=np.sqrt(x**2+y**2)
		mask=np.where(np.absolute(mid-a.z)<0.1*R)
		B=a.bfield[:,2] * code_units.B_cu *1e6	
		plt.figure()
		plt.scatter(np.log10(rs[mask]),np.log10(B[mask]),s=0.1)
		plt.xlabel('log10(r)'),plt.ylabel('log10(B_z)')


def histready(x,y,weight):
	hist_weighted, xb, yb = np.histogram2d(y, x, weights=weight, bins=(500, 500))
	hist_numbers, xb, yb = np.histogram2d(y, x, bins=(500, 500))
	hist_final = hist_weighted / hist_numbers
	hist_final = np.ma.masked_where(hist_numbers < 1, hist_final)
	ip = np.where(hist_numbers > 0)
	max_image = np.max(hist_final[ip])
	min_image = np.min(hist_final[ip])
	if ( (max_image/min_image) > 50 ):
		hist_final = np.nan_to_num(np.log10(hist_final))
	return hist_final,xb,yb

def sixframes(dirname,snapshot):
	a=arepo_utils.aread(dirname + snapshot[0])
	b=arepo_utils.aread(dirname + snapshot[1])
	c=arepo_utils.aread(dirname + snapshot[2])
	d=arepo_utils.aread(dirname + snapshot[3])
	e=arepo_utils.aread(dirname + snapshot[4])
	f=arepo_utils.aread(dirname + snapshot[5])
	fig,axs=plt.subplots(2,3,sharey=True,sharex=True)

	hist_final,xb,yb = histready(a.x,a.z,a.rho)
	axs[0,0].imshow(hist_final, aspect='auto', cmap='plasma', origin='lower', extent=[yb[0],yb[-1],xb[0],xb[-1]])

	hist_final,xb,yb = histready(b.x,b.z,b.rho)
	axs[0,1].imshow(hist_final, aspect='auto', cmap='plasma', origin='lower', extent=[yb[0],yb[-1],xb[0],xb[-1]])

	hist_final,xb,yb = histready(c.x,c.z,c.rho)
	axs[0,2].imshow(hist_final, aspect='auto', cmap='plasma', origin='lower', extent=[yb[0],yb[-1],xb[0],xb[-1]])

	hist_final,xb,yb = histready(d.x,d.z,d.rho)
	axs[1,0].imshow(hist_final, aspect='auto', cmap='plasma', origin='lower', extent=[yb[0],yb[-1],xb[0],xb[-1]])

	hist_final,xb,yb = histready(e.x,e.z,e.rho)
	axs[1,1].imshow(hist_final, aspect='auto', cmap='plasma', origin='lower', extent=[yb[0],yb[-1],xb[0],xb[-1]])

	hist_final,xb,yb = histready(f.x,f.z,f.rho)
	axs[1,2].imshow(hist_final, aspect='auto', cmap='plasma', origin='lower', extent=[yb[0],yb[-1],xb[0],xb[-1]])
