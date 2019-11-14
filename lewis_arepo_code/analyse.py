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

def centered(a,center,boxsize):
	'''returns midpoint and 3d distance to center
	center='mid' or 'rho', 'rho' gives density peak position'''
	if center=='mid':
		mid=boxsize/2
		midx,midy,midz=mid,mid,mid
	if center=='rho':
		mid=np.where(a.rho==a.rho.max())
		midx=a.x[mid]
		midy=a.y[mid]
		midz=a.z[mid]
	rs=np.sqrt((midx-a.x)**2+(midy-a.y)**2+(midz-a.z)**2)
	return midx,midy,midy,rs
	
def cylinder_velocity(a,center,boxsize):
	'''returns radial and rotational velocities (cylindrical)
	center='mid' or 'rho'''
	midx,midy,midz,rs=centered(a,center,boxsize)
	x,y=(a.x-midx),(a.y-midy)
	r=np.sqrt(x**2+y**2)
	vr=(x*a.vx+y*a.vy)/r
	vr=vr*code_units.v_cu  /1e5
	vtheta=(x*a.vy-y*a.vx)/(x**2+y**2)
	vtheta=vtheta*code_units.v_cu /1e5
	return vr,vtheta

def slice(a,boxsize,center,thickness):
	'''returns mask of elements within an x/y slice with z thickness,
	alsoe returns 3d distances from centre.
	centre=='mid' or 'rho', 'rho' makes centre where density peaks.'''
	midx,midy,midz,rs=centered(a,center,boxsize)
	h=np.sqrt((midz-a.z)**2)
	mask=np.where(h<=thickness)
	return mask,rs[mask]

def radial_profile(a,boxsize,center,thickness,weight_type):
	'''returns radius and radial profile within z slice,
	centre=='mid' or 'rho', 'rho' gives peak density position
	weight_type='rho'/'vr'/'vtheta' '''	
	mask,rs=slice(a,boxsize,center,thickness)
	x,y,z=a.x[mask],a.y[mask],a.z[mask]
	
	if weight_type=='rho':
		weight=a.rho[mask]
	if weight_type=='vr':
		weight=cylinder_velocity(a,center,boxsize)[0][mask]
	if weight_type=='vtheta':
		weight=cylinder_velocity(a,center,boxsize)[1][mask]
	return rs,weight

def avline(r,w):
	d=binned_statistic(r,w,bins=50)
	return d
	
def timeplot(dirname,snaps,boxsize,weight_type,ax,log):
	if weight_type=='rho':
		unit=code_units.rho_cu
		tag='log10(rho/gcm^-3)'
	if weight_type=='vr':
		unit=code_units.v_cu /1e5
		tag='kms^-1'
	if weight_type=='vtheta':
		unit=code_units.v_cu /1e5
		tag='kms^-1'
	for i in range(len(snaps)):
		a=arepo_utils.aread(dirname+snaps[i])
		rs,w=radial_profile(a,boxsize,'mid',boxsize/4*1e-3,weight_type)
		if log=='yes':
			d=avline(np.log10(rs*code_units.d_cu/ap.pc.cgs.value),np.log10(w*unit))
			ax.plot(d[1][:-1],d[0])
			ax.set_ylabel('%s.'%tag)
			ax.set_xlabel('log10(r/pc')
		else:
			d=avline(np.log10(rs*code_units.d_cu/ap.pc.cgs.value),w*unit)
			ax.plot(d[1][:-1],d[0])
			ax.set_ylabel('%s.'%tag)
			ax.set_xlabel('log10(r/pc)')
		
def multiplot(dirname,snaps,boxsize):
	fig,ax=plt.subplots(3,1)
	timeplot(dirname,snaps,boxsize,'rho',ax[0],'yes')
	timeplot(dirname,snaps,boxsize,'vr',ax[1],'no')
	timeplot(dirname,snaps,boxsize,'vtheta',ax[2],'no')
	return 	fig



def ratioB(a,boxsize):
	mid=boxsize/2
	bx=a.bfield[:,0]
	by=a.bfield[:,1]
	bz=a.bfield[:,2]
	B=np.sqrt(bx**2+by**2+bz**2)
	x=a.x-mid
	y=a.y-mid
	z=a.z-mid
	Btoroid=np.absolute(x*by-y*bx) #y*bz-z*by+z*bx-x*bz+xby-ybx
	Bpoloid=np.sqrt(B**2-Btoroid**2)
	ratio=Btoroid/Bpoloid
	radial=x*bx+y*by+z*bz
	print('max radial: '+str(radial.max()))
	return ratio



def histready(x,y,weight,force_lin):
	x=x*code_units.d_cu
	y=y*code_units.d_cu
	
		
	
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
	if force_lin != 'yes':
		if ( (max_image/min_image) > 50 ):
			print('log scale')
			hist_final = np.nan_to_num(np.log10(hist_final))
	return hist_final,xb,yb

def crop(a,zoomzone,boxsize):
	mid=boxsize/2
	mask1=np.where(np.absolute(mid-a.x)<zoomzone)
	mask2=np.where(np.absolute(mid-a.y)<zoomzone)
	mask3=np.where(np.absolute(mid-a.z)<zoomzone)
	MASK=np.intersect1d(mask1,mask2)
	MASK=np.intersect1d(MASK,mask3)
	return MASK	

def ploter(a,x,y,weight,zoomzone):
	mid=np.where(a.rho==a.rho.max())
	mask=np.where(np.absolute(a.x-a.x[mid])<zoomzone)
	mask1=np.where(np.absolute(a.y-a.y[mid])<zoomzone)
	mask2=np.where(np.absolute(a.z-a.z[mid])<zoomzone)
	MASK=np.intersect1d(mask,mask1)
	MASK=np.intersect1d(MASK,mask2)
	hist_final,xb,yb = histready(x[MASK],y[MASK],weight[MASK])		
	plt.imshow(hist_final, aspect='auto', cmap='plasma', origin='lower', extent=[yb[0],yb[-1],xb[0],xb[-1]])	




def plot6(dirname,snaps,weight_type,zoomzone,boxsize,force_lin):
	'''takes directory and array of 6 snapshot names, weight_type='rho' or 'bratio',zoomzone is maximum 
	radius from center in code units, boxsize in code units, force_lin='yes' if you want to block the 
	log scaling of plots.'''
	fig,axs=plt.subplots(2,3,sharex=True,sharey=True)
	for i in range (len(snaps)):
		a=arepo_utils.aread(dirname+snaps[i])
		if weight_type=='rho':
			weight=a.rho
			weightunit=code_units.rho_cu
		if weight_type=='bratio':
			weight=ratioB(a,boxsize)
			weightunit=1
		MASK=crop(a,zoomzone,boxsize)
		hist_final,xb,yb=histready(a.x[MASK]*code_units.d_cu,a.z[MASK]*code_units.d_cu,weight[MASK]*weightunit,force_lin)
		if i==0:
			im=axs[0,0].imshow(hist_final,aspect='auto', cmap='plasma', origin='lower', extent=[yb[0],yb[-1],xb[0],xb[-1]])
			clim=im.properties()['clim']
			axs[0,0].set_xticks([])
			axs[0,0].set_yticks([])
		if i==1:
			axs[0,1].imshow(hist_final,clim=clim,aspect='auto', cmap='plasma', origin='lower', extent=[yb[0],yb[-1],xb[0],xb[-1]])
			axs[0,1].set_xticks([])
			axs[0,1].set_yticks([])
		if i==2:
			axs[0,2].imshow(hist_final,clim=clim,aspect='auto', cmap='plasma', origin='lower', extent=[yb[0],yb[-1],xb[0],xb[-1]])
			axs[0,2].set_xticks([])
			axs[0,2].set_yticks([])             
		if i==3:
			axs[1,0].imshow(hist_final,clim=clim,aspect='auto', cmap='plasma', origin='lower', extent=[yb[0],yb[-1],xb[0],xb[-1]])
			axs[1,0].set_xticks([])
			axs[1,0].set_yticks([])                
		if i==4:
			axs[1,1].imshow(hist_final,clim=clim,aspect='auto', cmap='plasma', origin='lower', extent=[yb[0],yb[-1],xb[0],xb[-1]])
			axs[1,1].set_xticks([])
			axs[1,1].set_yticks([])                
		if i==5:
			axs[1,2].imshow(hist_final,clim=clim,aspect='auto', cmap='plasma', origin='lower', extent=[yb[0],yb[-1],xb[0],xb[-1]])
			axs[1,2].set_xticks([])
			axs[1,2].set_yticks([])
		fig.subplots_adjust(0.1,0.1,0.9,0.9,0,0)
		fig.colorbar(im,ax=axs.ravel().tolist(), shrink=1,pad=0)



