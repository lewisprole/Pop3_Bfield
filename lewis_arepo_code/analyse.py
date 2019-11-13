import arepo_utils 
import numpy as np
import glob
import matplotlib.pyplot as plt
import code_units 
import astropy.constants as ap
from scipy.stats import binned_statistic
import io
import sys
from scipy.stats import binned_statistic_2d

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
		DIF=1e100
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


def histready(x,y,weight,weight_type):
	x=x*code_units.d_cu
	y=y*code_units.d_cu
	if weight_type=='rho':
		weight=weight*code_units.rho_cu
		print(weight[0])
	hist_weighted, xb, yb = np.histogram2d(y, x, weights=weight, bins=(400, 400))
	
	hist_numbers, xb, yb = np.histogram2d(y, x, bins=(400, 400))
	
	#hist_mean=binned_statistic_2d(y,x,weight,bins=(400,400))[0]
	#hist_mean=hist_weighted/hist_numbers
	#MEAN=np.mean(weight)
	
	hist_final = hist_weighted / hist_numbers
	hist_final = np.ma.masked_where(hist_numbers < 1, hist_final)
	ip = np.where(hist_numbers > 0)
	max_image = np.max(hist_final[ip])
	min_image = np.min(hist_final[ip])
	if ( (max_image/min_image) > 50 ):
		hist_final = np.nan_to_num(np.log10(hist_final))
	return hist_final,xb,yb

def crop(a,zoomzone,provide_mid):
	if isinstance(provide_mid,str):
		mid=np.where(a.rho==a.rho.max())
		mask1=np.where(np.absolute(a.x[mid]-a.x)<zoomzone)
		mask2=np.where(np.absolute(a.y[mid]-a.y)<zoomzone)
		mask3=np.where(np.absolute(a.z[mid]-a.z)<zoomzone)
		MASK=np.intersect1d(mask1,mask2)
		MASK=np.intersect1d(MASK,mask3)
		return  MASK,(a.x[mid],a.y[mid],a.z[mid])
	else:
		mid=provide_mid
		mask1=np.where(np.absolute(mid[0]-a.x)<zoomzone)
		mask2=np.where(np.absolute(mid[1]-a.y)<zoomzone)
		mask3=np.where(np.absolute(mid[2]-a.z)<zoomzone)
		MASK=np.intersect1d(mask1,mask2)
		MASK=np.intersect1d(MASK,mask3)
		return MASK,mid	

def ploter(a,x,y,weight,zoomzone):
	mid=np.where(a.rho==a.rho.max())
	mask=np.where(np.absolute(a.x-a.x[mid])<zoomzone)
	mask1=np.where(np.absolute(a.y-a.y[mid])<zoomzone)
	mask2=np.where(np.absolute(a.z-a.z[mid])<zoomzone)
	MASK=np.intersect1d(mask,mask1)
	MASK=np.intersect1d(MASK,mask2)
	hist_final,xb,yb = histready(x[MASK],y[MASK],weight[MASK])		
	plt.imshow(hist_final, aspect='auto', cmap='plasma', origin='lower', extent=[yb[0],yb[-1],xb[0],xb[-1]])	

def sixframes(dirname,snapshot,times,zoomzone):
	a=arepo_utils.aread(dirname + snapshot[0])
	b=arepo_utils.aread(dirname + snapshot[1])
	c=arepo_utils.aread(dirname + snapshot[2])
	d=arepo_utils.aread(dirname + snapshot[3])
	e=arepo_utils.aread(dirname + snapshot[4])
	f=arepo_utils.aread(dirname + snapshot[5])
	fig,axs=plt.subplots(2,3,sharey=True,sharex=True)
	

	MASK,mid=crop(a,zoomzone,'no')
	hist_final,xb,yb = histready(a.x[MASK],a.z[MASK],a.rho[MASK],'rho')
	im=axs[0,0].imshow(hist_final, aspect='auto', cmap='plasma', origin='lower', extent=[yb[0],yb[-1],xb[0],xb[-1]],label=('%.f t_ff'%times[0]))
	clim=im.properties()['clim']
	axs[0,0].set_xticks([])
	axs[0,0].set_yticks([])
	axs[0,0].legend('%.f t_ff'%times[0])

	MASK,mid=crop(b,zoomzone,mid)
	hist_final,xb,yb = histready(b.x[MASK],b.z[MASK],b.rho[MASK],'rho')
	axs[0,1].imshow(hist_final, clim=clim,aspect='auto', cmap='plasma', origin='lower', extent=[yb[0],yb[-1],xb[0],xb[-1]],label=times[1])
	axs[0,1].set_xticks([])
	axs[0,1].set_yticks([])	
	axs[0,1].legend()

	MASK,mid=crop(c,zoomzone,mid)
	hist_final,xb,yb = histready(c.x[MASK],c.z[MASK],c.rho[MASK],'rho')
	im=axs[0,2].imshow(hist_final, clim=clim,aspect='auto', cmap='plasma', origin='lower', extent=[yb[0],yb[-1],xb[0],xb[-1]],label=times[2])
	axs[0,2].set_xticks([])
	axs[0,2].set_yticks([])
	axs[0,2].legend()

	MASK,mid=crop(d,zoomzone,mid)
	hist_final,xb,yb = histready(d.x[MASK],d.z[MASK],d.rho[MASK],'rho')
	im=axs[1,0].imshow(hist_final, clim=clim,aspect='auto', cmap='plasma', origin='lower', extent=[yb[0],yb[-1],xb[0],xb[-1]],label=times[3])
	axs[1,0].set_xticks([])
	axs[1,0].set_yticks([])
	axs[1,0].legend()

	MASK,mid=crop(e,zoomzone,mid)
	hist_final,xb,yb = histready(e.x[MASK],e.z[MASK],e.rho[MASK],'rho')
	im=axs[1,1].imshow(hist_final, clim=clim,aspect='auto', cmap='plasma', origin='lower', extent=[yb[0],yb[-1],xb[0],xb[-1]],label=times[4])
	axs[1,1].set_xticks([])	
	axs[1,1].set_yticks([])
	axs[1,1].legend()

	MASK,mid=crop(f,zoomzone,mid)
	hist_final,xb,yb = histready(f.x[MASK],f.z[MASK],f.rho[MASK],'rho')
	im=axs[1,2].imshow(hist_final, clim=clim,aspect='auto', cmap='plasma', origin='lower', extent=[yb[0],yb[-1],xb[0],xb[-1]],label=times[5])
	axs[1,2].set_xticks([])
	axs[1,2].set_yticks([])
	axs[1,2].legend()	

	fig.subplots_adjust(0.1,0.1,0.9,0.9,0,0)
	
	fig.colorbar(im,ax=axs.ravel().tolist(), shrink=1,pad=0)
