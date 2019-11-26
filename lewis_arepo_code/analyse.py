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





'''//////////////general snapshot query functions///////////////'''

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
		





'''//////////////functions for Hannebelle et al. 2008////////////////'''	


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
	vtheta=((x*a.vy-y*a.vx)/(x**2+y**2))*np.sqrt(x**2+y**2)
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
	if weight_type=='Bz':
		weight=a.bfield[:,2][mask]
	return rs,weight

def avline(r,w,bins):
	mask=np.isnan(w)
	d=binned_statistic(r[~mask],w[~mask],bins=bins)
	return d

def ratioB(a,boxsize,slice_option):
	if slice_option=='yes':
		mask,rs=slice(a,boxsize,'rho',boxsize/4*1e-3)
	else:
		mask=np.where(a.x==a.x) #everywhere 
	midx,midy,midz,rs=centered(a,'rho',boxsize)
	bx=a.bfield[:,0][mask]
	by=a.bfield[:,1][mask]
	bz=a.bfield[:,2][mask]
	B=np.sqrt(bx**2+by**2+bz**2)

	x=a.x[mask]-midx
	y=a.y[mask]-midy
	z=a.z[mask]-midz
	#Btoroid=x*by-y*bx #y*bz-z*by+z*bx-x*bz+xby-ybx
	#Bpoloid=bz
	#ratio=np.absolute(Btoroid/Bpoloid)
	ratio=1/np.sqrt(1+(x/y)**2) *(bx-x/y*by) / bz
	return np.absolute(ratio)
	
def timeplot(dirname,snaps,boxsize,weight_type,ax,log,bins):
	if weight_type=='rho':
		unit=code_units.rho_cu
		tag='log10(rho/gcm^-3)'
		y1,y2=-18.5,-11.5
	if weight_type=='vr':
		unit=code_units.v_cu /1e5
		tag='v_r/kms^-1'
		y1,y2=-1.4,0
	if weight_type=='vtheta':
		unit=code_units.v_cu /1e5
		tag='v_theta/kms^-1'
		y1,y2=0,2
	if weight_type=='Bz':
		unit=code_units.B_cu*1e6
		tag='Bz/uG'
		y1,y2=1,6
	for i in range(len(snaps)):
		a=arepo_utils.aread(dirname+snaps[i])
		rs,w=radial_profile(a,boxsize,'rho',boxsize/4*1e-3,weight_type)
		if log=='yes':
			d=avline(np.log10(rs[np.where(rs>0)]*code_units.d_cu/ap.pc.cgs.value),np.log10(w[np.where(rs>0)]*unit),bins)
			ax.plot(d[1][:-1],d[0])
			ax.set_ylabel('%s.'%tag)
			ax.set_xlabel('log10(r/pc)')
			ax.set_xlim(-4.5,-1.5)
			ax.set_ylim(y1,y2)
		else:
			d=avline(np.log10(rs*code_units.d_cu/ap.pc.cgs.value),w*unit,bins)
			ax.plot(d[1][:-1],d[0])
			ax.set_ylabel('%s.'%tag)
			ax.set_xlabel('log10(r/pc)')
			ax.set_xlim(-4.5,-1.5)
			ax.set_ylim(y1,y2)

def multiplot(dirname,snaps,boxsize,mu,B):
	fig,ax=plt.subplots(4,1)
	ax[0].set_title('mu=%i.,B=%f.'%(mu,B))
	timeplot(dirname,snaps,boxsize,'rho',ax[0],'yes',45)
	timeplot(dirname,snaps,boxsize,'vr',ax[1],'no',45)
	timeplot(dirname,snaps,boxsize,'vtheta',ax[2],'no',45)
	timeplot(dirname,snaps,boxsize,'Bz',ax[3],'yes',30)
	return 	fig









'''//////////functions for Florian et al. 2011//////////////'''


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
	hist_final,xb,yb = histready(x[MASK],y[MASK],weight[MASK],'no')		
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
			weight=np.log10(ratioB(a,boxsize,'no'))
			weightunit=1
		if weight_type=='v':
			weight=np.sqrt(a.vx**2+a.vy**2+a.vz**2)
			weightunit=code_units.v_cu
		MASK=crop(a,zoomzone,boxsize)
		#mask=np.where(a.y<boxsize/2)
		#MASK=np.intersect1d(MASK,mask)
		hist_final,xb,yb=histready(a.x[MASK]*code_units.d_cu,a.z[MASK]*code_units.d_cu,weight[MASK]*weightunit,force_lin)
		ax=axs.ravel()
		if i==0:
			im=ax[0].imshow(hist_final,aspect='auto', cmap='plasma', origin='lower', extent=[yb[0],yb[-1],xb[0],xb[-1]])
			clim=im.properties()['clim']
			ax[0].set_xticks([])
			ax[0].set_yticks([])
		else:
			ax[i].imshow(hist_final,clim=clim,aspect='auto', cmap='plasma', origin='lower', extent=[yb[0],yb[-1],xb[0],xb[-1]])
			ax[i].set_xticks([])
			ax[i].set_yticks([])
		
		fig.subplots_adjust(0.1,0.1,0.9,0.9,0,0)
		fig.colorbar(im,ax=axs.ravel().tolist(), shrink=1,pad=0)


def sliceplot(dirname,names,weight_type):
	if weight_type=='rho':
		tag='log10(rho/gcm^-3)'
		unit=code_units.rho_cu
	if weight_type=='bratio':
		unit=code_units.B_cu
		tag='log10(Bz/G)'
	fig,axs=plt.subplots(2,3,sharey=True,sharex=True)
	axs=axs.ravel()
	for i in range(len(axs)):
		A=np.fromfile(dirname+names[i], dtype='int32', sep="")[2:].reshape([1000, 1000])
		if i==0:
			im=axs[i].imshow(np.log10(A.T*unit),aspect='auto', cmap='plasma', origin='lower')
			clim=im.properties()['clim']
			axs[i].set_xticks([])
			axs[i].set_yticks([])
			axs[i].set_xlim(0,1000)
			axs[i].set_ylim(0,1000)
		else:
			axs[i].imshow(np.log10(A.T*unit),clim=clim,aspect='auto', cmap='plasma', origin='lower')
			axs[i].set_xticks([])
			axs[i].set_yticks([])
			axs[i].set_xlim(0,1000)
			axs[i].set_ylim(0,1000)
	fig.subplots_adjust(0.1,0.1,0.9,0.9,0,0)
	cbar=fig.colorbar(im,ax=axs.tolist(), shrink=1,pad=0)
	cbar.ax.set_ylabel(tag, rotation=270,labelpad=25)

def ratioB_cube(name,size_cu,pixels):
	
	x=np.linspace(0,size_cu,pixels) 
	y,x=np.meshgrid(x,x)
	x=x-size_cu/2
	y=y-size_cu/2
	
	A=np.fromfile(name,dtype=np.int32)[3:]
	Bx=A[0::3].reshape(pixels,pixels,pixels)
	By=A[1::3].reshape(pixels,pixels,pixels)
	Bz=A[2::3].reshape(pixels,pixels,pixels)
	bratio=np.zeros((pixels,pixels,pixels))
	for i in range (pixels):
		bx=Bx[:,:,i]
		by=By[:,:,i]
		bz=Bz[:,:,i]
		bratio[:,:,i]=1/np.sqrt(1+(x/y)**2) * (bx-(x/y)*by) /bz
	#Brot=1/np.sqrt(1+(x/y)**2) * (Bx-(x/y)*By)
	#bratio=Brot/Bz
	return bratio

def ratioB_plot(dirname,names,weight_type,pixels):
	fig,axs=plt.subplots(2,3,sharey=True,sharex=True)
	axs=axs.ravel()
	for i in range (len(axs)):
		if weight_type=='bratio':
			bratio=np.absolute(ratioB_cube(dirname+names[5-i],0.12,pixels)[:,350:650,:])
			bratio=np.sum(bratio,1)
		if weight_type=='rho':
			print(i)
			bratio=np.fromfile(dirname+names[5-i],dtype=np.int32)[3:].reshape(pixels,pixels,pixels)
			bratio=np.sum(bratio[:,480:520,:],1)
		if i==0:
			im=axs[5-i].imshow(np.log10(np.rot90(bratio)))
			clim=im.properties()['clim']
		else:	
			axs[5-i].imshow(np.log10(np.rot90(bratio)),clim=clim)
		axs[5-i].set_xticks([])
		axs[5-i].set_yticks([])
		axs[5-i].set_xlim(0,1000)
		axs[5-i].set_ylim(0,1000)
	fig.subplots_adjust(0.1,0.1,0.9,0.9,0,0)
	cbar=fig.colorbar(im,ax=axs.tolist(), shrink=1,pad=0)
	cbar.ax.set_ylabel('log10(Brot/Bz)', rotation=270,labelpad=25)
	



def spacial_average(name,size_cu,pixels):
	'''look down on xy plane and average the B_rot/B_pol ratio for each 2d R'''
	A=np.fromfile(name,dtype=np.int32)
	bratio=np.absolute(ratioB_cube(name,0.12,pixels)[:,:,400:600])
	bratio=np.sum(bratio,2).flatten()
	x=np.linspace(-size_cu/2,size_cu/2,pixels)
	y,x=np.meshgrid(x,x)
	x=x.flatten()
	y=y.flatten()
	rs=np.sqrt(y**2+x**2)
	


	#a=arepo_utils.aread(snap)
	#midx,midy,midz,rs=centered(a,'rho',boxsize)	
	#bratio=ratioB(a,boxsize,'no')
	#rs=np.sqrt((a.x-midx)**2+(a.y-midy)**2)*code_units.d_cu*1e-17
	#mask,junk=slice(a,boxsize,'rho',boxsize/4*1e-3)
	#rs,bratio=rs,bratio#rs[mask],bratio[mask]
	mask=np.where(rs<0.05)
	d=avline(rs[mask],bratio[mask],20)
	return d[1][:-1],d[0]
	
def spacial_multiplot(dirname,snaps,im_size_cu,pixels):
	fig,ax=plt.subplots(1)
	fig.set_label('r/ 10^17 cm')
	fig.set_label('B_rot/B_pol')
	for i in range(len(snaps)):
		x,y=spacial_average(dirname+snaps[i],im_size_cu,pixels)
		ax.plot(x,y,label=snaps[i][-3:])
	ax.set_xlim(0,0.05)
	ax.legend()


