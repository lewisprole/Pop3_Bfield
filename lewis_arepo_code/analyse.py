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
from  mpl_toolkits.axes_grid1 import ImageGrid



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

'''/////////plotting functions using 2d histogram methods from cell positions//////////''' 

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

'''//////////plotting functions for arepo grid projection files/////////////'''

def read_cube(name):
	'''returns 1 cube[x,y,z] for a scaler vaiable'''
	A=np.fromfile(name,dtype=np.int32)
	n=A[0]
	N=n*n
	A=A[2:]
	cube=np.zeros((n,n,n))
	for i in range(n):
		for j in range(n):
			line=A[i*N+j*n:i*N+j*n+n]
			cube[i,j,:]=line
	return cube


def read_3cube(name):
	'''returns 3 cubes giving x,y and z components, individual cubes have positional
	axis cube[x,y,z]'''
	A=np.fromfile(name,dtype=np.int32)
	n=A[0]
	N=n*n
	A=A[3:]
	Ax=A[0::3]
	Ay=A[1::3]
	Az=A[2::3]
	cubex=np.zeros((n,n,n))
	cubey=np.zeros((n,n,n))
	cubez=np.zeros((n,n,n))
	for i in range(n):
		for j in range(n):
			linex=Ax[i*N+j*n:i*N+j*n+n]
			cubex[i,j,:]=linex
			liney=Ay[i*N+j*n:i*N+j*n+n]
			cubey[i,j,:]=liney
			linez=Az[i*N+j*n:i*N+j*n+n]
			cubez[i,j,:]=linez
	return cubex,cubey,cubez


def reduce_velocity_resolution(vx,vy,reduc):
	'''takes 2 2D arrays for x and y velocities, and a reduction factor
	returns 2 2D arrays with side length reduced by N/reduc by averaging velocities
	used for quiver plots'''
	
	n=int(len(vx[:,0])/reduc)
	velx,vely=np.zeros((n,n)),np.zeros((n,n))
	for i in range(n):
		for j in range(n):
			valx=np.mean(vx[reduc*i:reduc*i+reduc,reduc*j:reduc*j+reduc])
			velx[i,j]=valx
			valy=np.mean(vy[reduc*i:reduc*i+reduc,reduc*j:reduc*j+reduc])
			vely[i,j]=valy
	return velx,vely


def arrowplot(ax,vx,vz,reduc):
	'''xz plane quiver plot, vx and vz should be 2D arrays''' 
	
	pixels=len(vx[0,:])
	vx,vz=reduce_velocity_resolution(vx,vz,reduc)
	x=np.linspace(0,pixels,int(pixels/reduc))
	z,x=np.meshgrid(x,x)

	z,x=np.rot90(z),np.rot90(x) #confusing rotation stuff 
	temp=vx
	vx=vz
	vz=-temp	

	ax.quiver(z,x,vz,vx,headwidth=10,minshaft=3,pivot='mid',color='w')
	ax.set_ylim(0,pixels)
	ax.set_xlim(0,pixels)

def add_arrows(axs,dirname,names,width,reduc):
	'''addition to ratioB_plot function'''
	for i in range(len(axs)):
		print(i)
		vx,vy,vz=read_3cube(dirname+names[i])	
		mid=int(len(vx[0,0,:])/2)
		vx=np.sum(vx[:,mid-width:mid+width,:],1)
		vz=np.sum(vz[:,mid-width:mid+width,:],1)	
		arrowplot(axs[i],vx,vz,reduc)

def sliceplot(dirname,names,weight_type,pixels):
	if weight_type=='rho':
		tag='log10(rho/gcm^-3)'
		unit=code_units.rho_cu
	if weight_type=='bratio':
		unit=code_units.B_cu
		tag='log10(Bz/G)'
	fig,axs=plt.subplots(2,3,sharey=True,sharex=True)
	axs=axs.ravel()
	for i in range(len(axs)):
		A=np.fromfile(dirname+names[i], dtype='int32', sep="")[2:].reshape([pixels,pixels])
		if i==0:
			im=axs[i].imshow(np.log10(A.T*unit),aspect='auto', cmap='plasma', origin='lower')
			clim=im.properties()['clim']
			axs[i].set_xticks([])
			axs[i].set_yticks([])
			axs[i].set_xlim(0,pixels)
			axs[i].set_ylim(0,pixels)
		else:
			axs[i].imshow(np.log10(A.T*unit),clim=clim,aspect='auto', cmap='plasma', origin='lower')
			axs[i].set_xticks([])
			axs[i].set_yticks([])
			axs[i].set_xlim(0,pixels)
			axs[i].set_ylim(0,pixels)
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
		bratio[:,:,i]=(x*by-y*bx)/np.sqrt(x**2+y**2) /bz   #1/np.sqrt(1+(x/y)**2) * (bx-(x/y)*by) /bz
	#Brot=1/np.sqrt(1+(x/y)**2) * (Bx-(x/y)*By)
	#bratio=Brot/Bz
	return bratio

def ratioB_plot(dirname,names,weight_type,pixels):
	fig=plt.figure(figsize=(4,4))
	grid=ImageGrid(fig,111,nrows_ncols=(2,3),axes_pad=0,cbar_mode='single')#,axs=plt.subplots(2,3)#sharey=True,sharex=True)
	#axs=axs.ravel()
	width_rho=int(pixels/40)
	width_bratio=int(pixels/10)
	mid=int(pixels/2)
	for i in range(len(grid)):
		print(i)
		if weight_type=='bratio':
			bratio=np.absolute(ratioB_cube(dirname+names[5-i],0.12,pixels)[:,mid-width_bratio:mid+width_bratio,:])
			bratio=np.sum(bratio,1)/(2*width_bratio)
			tag='log10(Brot/Bz)'
		if weight_type=='rho':
			bratio=read_cube(dirname+names[5-i])*code_units.rho_cu#(np.fromfile(dirname+names[5-i],dtype=np.int32)[3:].reshape(pixels,pixels,pixels) *code_units.rho_cu
			bratio=np.sum(bratio[:,mid-width_rho:mid+width_rho,:],1)/(2*width_rho)
			tag='log10(rho/gcm^3)'
		if i==0:
			
			im=grid[5-i].imshow(np.log10(np.rot90(bratio)),cmap='plasma')
			clim=im.properties()['clim']
		else:	
			im=grid[5-i].imshow(np.log10(np.rot90(bratio)),clim=clim,cmap='plasma')
		#axs[5-i].set_ylim([0,pixels])
		#axs[5-i].set_xlim([0,pixels])
		grid[5-i].set_yticks([])	
		grid[5-i].set_xticks([])
	#fig.subplots_adjust(left=0.1,right=0.9,bottom=0.1,top=0.9,wspace=0,hspace=-0.4)
	Grid=list((grid[0],grid[1],grid[2],grid[3],grid[4],grid[5]))
	cbar = grid.cbar_axes[0].colorbar(im)
	
	#cbar=fig.colorbar(im,ax=grid,shrink=1,pad=0)
	cbar.ax.set_ylabel(tag, rotation=270,labelpad=25)
	return fig,grid,im

def rho_arrow_plot(dirname,rho_names,vel_names,pixels):
	fig,grid,im=ratioB_plot(dirname,rho_names,'rho',pixels)
	add_arrows(grid,dirname,vel_names,int(pixels/40),35)	

def spacial_average(name,size_cu,pixels):
	'''look down on xy plane and average the B_rot/B_pol ratio for each 2d R'''
	
	bratio=np.absolute(ratioB_cube(name,0.12,pixels)[:,:,480:520])
	bratio=np.sum(bratio,2).flatten()
	bratio=bratio/40
	x=np.linspace(-size_cu/2,size_cu/2,pixels)
	y,x=np.meshgrid(x,x)
	x=x.flatten()
	y=y.flatten()
	rs=np.sqrt(y**2+x**2)*code_units.d_cu*1e-17
	


	#a=arepo_utils.aread(snap)
	#midx,midy,midz,rs=centered(a,'rho',boxsize)	
	#bratio=ratioB(a,boxsize,'no')
	#rs=np.sqrt((a.x-midx)**2+(a.y-midy)**2)*code_units.d_cu*1e-17
	#mask,junk=slice(a,boxsize,'rho',boxsize/4*1e-3)
	#rs,bratio=rs,bratio#rs[mask],bratio[mask]
	mask=np.where(rs<=0.05)
	rs=rs[mask]
	bratio=bratio[mask]
	d=avline(rs,bratio,50)
	return d[1][:-1],d[0]
	
def spacial_multiplot(dirname,snaps,im_size_cu,pixels):
	fig,ax=plt.subplots(1)
	fig.set_label('r/ 10^17 cm')
	fig.set_label('B_rot/B_pol')
	for i in range(len(snaps)):
		x,y=spacial_average(dirname+snaps[i],im_size_cu,pixels)
		ax.plot(x,y,label=snaps[i][-3:])
	ax.semilogy()
	ax.set_xlim(0,0.05)
	ax.legend()


