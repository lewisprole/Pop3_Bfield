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
from scipy.interpolate import interp1d


'''//////////////general snapshot query functions///////////////'''

def sinkcheck(dirname):
	'''looks at all snapshots and tells you when the first sink particle is formed'''
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
	'''cycles through snapshots to find closest file to target time (seconds)'''
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
	'''returns mask of elements within an xy slice of dz thickness,
	also returns 3d distances from centre.
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
	'''average line through data'''
	mask=np.isnan(w)
	d=binned_statistic(r[~mask],w[~mask],bins=bins)
	return d

def ratioB(a,boxsize,slice_option,dz):
	'''returns ratio of rotational Bfield in xy pane to z Bfield
	slice_option=='yes' for a thin slice dz around z=mid'''
	if slice_option=='yes':
		mask,rs=slice(a,boxsize,'rho',dz)
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
	ratio=1/np.sqrt(1+(x/y)**2) *(bx-x/y*by) / bz
	RS=np.sqrt(x**2+y**2)

	return np.absolute(ratio),RS
	
def timeplot(dirname,snaps,boxsize,weight_type,ax,log,bins):
	'''radial profile of variable within thin slice in z direction, e.g. 'rho' for density'''
	t_ff=3e4
	if weight_type=='rho': #density 
		unit=code_units.rho_cu
		tag='log10(rho/gcm^-3)'
		y1,y2=-18.5,-11.5
	if weight_type=='vr': #radial velocity 
		unit=code_units.v_cu /1e5
		tag='v_r/kms^-1'
		y1,y2=-1.4,0
	if weight_type=='vtheta': #rotational velocity 
		unit=code_units.v_cu /1e5
		tag='v_theta/kms^-1'
		y1,y2=0,2
	if weight_type=='Bz': #z component of B field
		unit=code_units.B_cu*1e6
		tag='Bz/uG'
		y1,y2=1,6
	for i in range(len(snaps)):
		a=arepo_utils.aread(dirname+snaps[i])
		t=a.time*code_units.t_cu/(60*60*24*365) / t_ff
		rs,w=radial_profile(a,boxsize,'rho',boxsize/4*1e-3,weight_type)
		if log=='yes':
			d=avline(np.log10(rs[np.where(rs>0)]*code_units.d_cu/ap.pc.cgs.value),np.log10(w[np.where(rs>0)]*unit),bins)
			ax.plot(d[1][:-1],d[0],label='%.3f t/t_ff'%t)
			ax.set_ylabel('%s.'%tag)
			ax.set_xlabel('log10(r/pc)')
			ax.set_xlim(-4.5,-1.5)
			ax.set_ylim(y1,y2)
		else:
			d=avline(np.log10(rs*code_units.d_cu/ap.pc.cgs.value),w*unit,bins)
			ax.plot(d[1][:-1],d[0],label='%.3f t/t_ff'%t)
			ax.set_ylabel('%s.'%tag)
			ax.set_xlabel('log10(r/pc)')
			ax.set_xlim(-4.5,-1.5)
			ax.set_ylim(y1,y2)
	return ax	


def multiplot(dirname,snaps,boxsize,mu,B):
	'''6 pannel radial profile graphs'''
	fig,ax=plt.subplots(4,1)
	ax[0].set_title('mu=%i.,B=%f.'%(mu,B))
	ax[0]=timeplot(dirname,snaps,boxsize,'rho',ax[0],'yes',45)
	ax[0].legend(loc='upper right')
	timeplot(dirname,snaps,boxsize,'vr',ax[1],'no',45)
	timeplot(dirname,snaps,boxsize,'vtheta',ax[2],'no',45)
	timeplot(dirname,snaps,boxsize,'Bz',ax[3],'yes',30)
	
	return 	fig


def xy_Bratio(dirname,names,bins):
	for i in range(len(names)):
		print(i)
		bratio=np.absolute(ratioB_cube(dirname+names[i],0.06,600)[:,:,290:310])
		bratio=np.sum(bratio,2)/20
		x=np.linspace(-0.06,0.06,600)
		y,x=np.meshgrid(x,x)
		rs=np.sqrt(x**2+y**2)
		
		
		


		#a=arepo_utils.aread(dirname+names[i])
		#t=a.time*code_units.t_cu/(60*60*24*265)/1e4
		#bratio,rs=ratioB(a,1.975,'yes',dz)
		mask=np.where(rs*code_units.d_cu/1e17<=0.054)
		d=avline(rs[mask],bratio[mask],bins)
		plt.plot(d[1][:-1]*code_units.d_cu/1e17,d[0])#,label='%.3f t_ff'%t)
		plt.xlim(0,0.05)
		plt.yscale('log')
	plt.legend(loc='upper right')
			






'''//////////functions for Florian et al. 2011//////////////'''

'''/////////plotting functions using 2d histogram methods from cell positions//////////''' 

def histready(x,y,weight,force_lin):
	'''2D histogram with weighting e.g. density ('rho')'''
	x=x*code_units.d_cu
	y=y*code_units.d_cu
	
		
	
	hist_weighted, xb, yb = np.histogram2d(y, x, weights=weight, bins=(400, 400))
	
	hist_numbers, xb, yb = np.histogram2d(y, x, bins=(400, 400))
	#hist_weighted,xb,yb=plt.hist2d(x,y,weights=weight,bins=[400,400])
	#hist_numbers,xb,yb=plt.hist2d(x,y,bins=[400,400])
	
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
	'''crops snapshot data to within a distance=zoomzone from the centre (in code units)'''
	mid=np.where(a.rho==a.rho.max())
	mask1=np.where(np.absolute(a.x[mid]-a.x)<zoomzone)
	mask2=np.where(np.absolute(a.y[mid]-a.y)<zoomzone)
	mask3=np.where(np.absolute(a.z[mid]-a.z)<zoomzone)
	MASK=np.intersect1d(mask1,mask2)
	MASK=np.intersect1d(MASK,mask3)
	return MASK	

def ploter(a,x,y,weight,zoomzone):
	'''plots cropped 2D histogram'''
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
			weight,j=np.log10(ratioB(a,boxsize,'no',0))
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





'''//////////Reading grid projection files/////////////'''

def read_cube(name):
	'''returns 1 cube[x,y,z] for a scaler vaiable'''
	A=np.fromfile(name,dtype=np.float32)
	a=np.fromfile(name,dtype=np.int32)
	n=a[0] #length of cube side #
	N=n*n #size of 2d face #
	A=A[2:] #cube data #
	cube=np.zeros((n,n,n)) #create cube space 

	for i in range(n):
		for j in range(n):
			line=A[i*N+j*n:i*N+j*n+n]
			cube[i,j,:]=line
	return cube


def read_3cube(name):
	'''returns 3 cubes giving x,y and z components, individual cubes have positional
	axis cube[x,y,z]'''
	A=np.fromfile(name,dtype=np.float32)
	a=np.fromfile(name,dtype=np.int32)
	n=a[0]
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



'''//////////Velocity arrow plots/////////////'''

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

#	n=vx.shape[0]
#	N=int(n/reduc)
#	velx,vely=vx[0::N,0::N],vy[0::N,0::N]
	
	return velx,vely


def arrowplot(ax,vx,vz,reduc,crop):
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

def add_arrows(axs,dirname,names,width,reduc,color,crop):
	'''addition to ratioB_plot function'''
	for i in range(len(axs)):
		print(i)
		vx,vy,vz=read_3cube(dirname+names[i])	
		mid=int(len(vx[0,0,:])/2)
		vx=vx[:,300,:]#np.sum(vx[:,mid-width:mid+width,:],1)/(2*width)
		vz=vz[:,300,:]#np.sum(vz[:,mid-width:mid+width,:],1)/(2*width)	
		arrowplot(axs[i],vx,vz,reduc,crop)

def arrows_version2(ax,dirname,vel_name,spacing,crop):
	vx,vy,vz=read_3cube(dirname+vel_name)
	pixels=len(vx[0,0,:])
	mid=int(pixels/2)
	width=int(pixels/40)
	vx=np.sum(vx[:,mid-width:mid+width,:],1)/(2*width)#[:,300,:]
	vz=np.sum(vz[:,mid-width:mid+width,:],1)/(2*width)#300,:]
	

	v=np.sqrt(vx**2+vz**2) #normalise to unit vector 
	vx,vz=vx/v,vz/v

	x=np.linspace(0,pixels,pixels)
	z,x=np.meshgrid(x,x)

	vx=vx[mid-crop:mid+crop,mid-crop:mid+crop] #crop the image for better arrow sizes
	vz=vz[mid-crop:mid+crop,mid-crop:mid+crop]
	x=x[mid-crop:mid+crop,mid-crop:mid+crop]
	z=z[mid-crop:mid+crop,mid-crop:mid+crop]

	x=x[0::spacing,0::spacing] #reduce number of arrows 
	z=z[0::spacing,0::spacing]
	vx=vx[0::spacing,0::spacing]
	vz=vz[0::spacing,0::spacing]

	ax.quiver(x[1:,1:],z[1:,1:],vx[1:,1:],vz[1:,1:],headwidth=5,minshaft=3,pivot='mid')#,scale=100)
	ax.set_ylim(mid-crop,mid+crop)
	ax.set_xlim(mid-crop,mid+crop)



'''//////////Plotting slice AREPO files/////////////'''

def sliceplot(dirname,names,weight_type,pixels):
	if weight_type=='rho':
		tag='log10(rho/gcm^-3)'
		unit=code_units.rho_cu
	if weight_type=='B':
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





'''//////////working with magnetic grid projection files/////////////'''

def ratioB_cube(name,zoomzone,pixels):
	'''ratio of toridal over poloidal B field'''	
	x=np.linspace(-zoomzone,zoomzone,pixels) 
	y,x=np.meshgrid(x,x)
	Bx,By,Bz=read_3cube(name)	
	bratio=np.zeros((pixels,pixels,pixels))
	for i in range (pixels):
		bx=Bx[:,:,i]
		by=By[:,:,i]
		bz=Bz[:,:,i]
		bratio[:,:,i]=(x*by-y*bx)/np.sqrt(x**2+y**2) /bz   

	return bratio

def B_vs_z(name,r,imradius_CU,pixels):
	'''plots arrays of average |B| vs z for a given r
	r is given in code units, if r is an array, each r is
	 plotted as a different line'''
	x=np.linspace(-imradius_CU,imradius_CU,pixels)
	y,x=np.meshgrid(x,x)
	rs=np.sqrt(y**2+x**2)
	pixel_size=2*imradius_CU/pixels
	dr=pixel_size*5
	Bx,By,Bz=read_3cube(name)
	B=np.sqrt(Bx**2+By**2+Bz**2)
	plt.figure()
	for j in range(len(r)):
		Bs=np.array([])
		zs=np.array([])
		for i in range(pixels):
			B_face=B[:,:,i]
			r_mask=np.where((rs<r[j]+dr) & (rs>r[j]-dr))
			B_value=(B_face[r_mask])
			N=len(B_value)
			B_value=sum(B_value)/N #average value in dr range 
			z=i*pixel_size-(pixels/2*pixel_size) #want the middle to be z=0
			Bs=np.append(Bs,B_value*code_units.B_cu*1e3) #B in mG 
			zs=np.append(zs,z*code_units.d_cu/ap.au.cgs.value) #z in AU
		plt.plot(Bs,zs,label='r=%.0f AU'%(r[j]*code_units.d_cu/ap.au.cgs.value)),plt.ylabel('z / AU'),plt.xlabel('B / mG')
	plt.legend()
	



'''//////////Plotting flattened grid projection files/////////////'''

def cube_plot(dirname,names,weight_type,pixels):
	'''set of 6 subplots from snapshot names, can choose weight_type as 'bratio' for 
	B_tor/B_pol ratio or 'rho' for a simple density display. Pixels is the grid projection
	side dimension. '''

	fig=plt.figure(figsize=(4,4))
	grid=ImageGrid(fig,111,nrows_ncols=(2,3),axes_pad=0,cbar_mode='single')
	width_rho=int(pixels/40)
	width_bratio=int(pixels/10)
	mid=int(pixels/2)
	for i in range(len(grid)):
		print(i)

		if weight_type=='bratio':
			bratio=np.absolute(ratioB_cube(dirname+names[5-i],0.06,pixels)[:,mid-width_bratio:mid+width_bratio,:])
			bratio=np.sum(bratio,1)/(2*width_bratio)
			tag='log10(Brot/Bz)'

		if weight_type=='rho':
			bratio=read_cube(dirname+names[5-i])*code_units.rho_cu
			bratio=np.sum(bratio[:,mid-width_rho:mid+width_rho,:],1)/(2*width_rho)
			tag='log10(rho/gcm^3)'
		if weight_type=='v':
			vx,vy,vz=read_3cube(dirname+names[5-i])
			bratio=np.sqrt(vx**2+vy**2+vz**2)*code_units.v_cu/1e5
			bratio=np.sum(bratio[:,mid-width_rho:mid+width_rho,:],1)/(2*width_rho)
			bratio=10**(bratio)
			tag='v/kms^-1'
		if i==0:	
			im=grid[5-i].imshow(np.log10(np.rot90(bratio)),cmap='plasma')
			clim=im.properties()['clim']
		else:	
			im=grid[5-i].imshow(np.log10(np.rot90(bratio)),clim=clim,cmap='plasma')
		
		grid[5-i].set_yticks([])	
		grid[5-i].set_xticks([])

	Grid=list((grid[0],grid[1],grid[2],grid[3],grid[4],grid[5]))
	cbar = grid.cbar_axes[0].colorbar(im)	
	cbar.ax.set_ylabel(tag, rotation=270,labelpad=25)
	return fig,grid,im


def plot3x3(dirname,snapnumbers,zoomzone,pixels):
	'''3 rho plots, 3 velocity plots, 3 bratio plots'''
	fig=plt.figure(figsize=(4,4))
	grid=ImageGrid(fig,111,nrows_ncols=(3,3),axes_pad=0,cbar_mode='edge')
	for i in range (len(snapnumbers)):
		rho=read_cube(dirname+'im'+snapnumbers[2-i]+'/density_grid_'+snapnumbers[2-i])[150:450,285:315,150:450] *code_units.rho_cu
		print('rho' + snapnumbers[2-i] + 'read')
		vx,vy,vz=read_3cube(dirname+'im'+snapnumbers[2-i]+'/velocity_grid_'+snapnumbers[2-i])
		v=np.sqrt(vx**2+vy**2+vz**2)
		v_xz=np.sqrt(vx**2+vz**2)
		vx,vz=(vx/v_xz)[150:450,300,150:450],(vz/v_xz)[150:450,300,150:450]
		v=v[150:450,285:315,150:450]*code_units.v_cu/1e5
		x=np.linspace(0,int(pixels/2)-1,int(pixels/2))
		z,x=np.meshgrid(x,x)

		print('v' + snapnumbers[2-i] + 'read')
		bratio=np.absolute(ratioB_cube(dirname+'im'+snapnumbers[2-i]+'/magnetic_grid_'+snapnumbers[2-i],zoomzone,pixels))[150:450,285:315,150:450]
		print('B' + snapnumbers[2-i] + 'read')

		if i==0: #setting up the colour scales 
			im_rho=grid[2-i].imshow(np.log10(np.rot90(np.sum(rho,1)/30)),cmap='plasma')
			clim_rho=im_rho.properties()['clim']
			cbar = grid.cbar_axes[0].colorbar(im_rho)
			cbar.ax.set_ylabel('rho/gcm^-3', rotation=270,labelpad=25)

			im_vel=grid[5-i].imshow(np.rot90(np.sum(v,1)/30),cmap='plasma')
			clim_vel=im_vel.properties()['clim']
			cbar = grid.cbar_axes[1].colorbar(im_vel)
			cbar.ax.set_ylabel('v/kms^-1', rotation=270,labelpad=25)
			grid[5-i].quiver(x[0::19,0::19][1:,1:],z[0::19,0::19][1:,1:],vx[0::19,0::19][1:,1:],vz[0::19,0::19][1:,1:],color='c',headwidth=4,pivot='mid')
			grid[5-i].set_ylim(0,300)


			im_b=grid[8-i].imshow(np.log10(np.rot90(np.sum(bratio,1)/30)),cmap='plasma')
			clim_b=im_b.properties()['clim']
			cbar = grid.cbar_axes[2].colorbar(im_b)
			cbar.ax.set_ylabel('B_tor/B_pol', rotation=270,labelpad=25)
		else:
			grid[2-i].imshow(np.log10(np.rot90(np.sum(rho,1)/30)),clim=clim_rho,cmap='plasma')
			grid[5-i].imshow(np.rot90(np.sum(v,1)/30),clim=clim_vel,cmap='plasma')
			grid[5-i].quiver(x[0::19,0::19][1:,1:],z[0::19,0::19][1:,1:],vx[0::19,0::19][1:,1:],vz[0::19,0::19][1:,1:],color='c',headwidth=4,pivot='mid')
			grid[5-i].set_ylim(0,300)
			grid[8-i].imshow(np.log10(np.rot90(np.sum(bratio,1)/30)),clim=clim_b,cmap='plasma')

		grid[2-i].set_yticks([])
		grid[2-i].set_xticks([])
		grid[5-i].set_yticks([])
		grid[5-i].set_xticks([])
		grid[8-i].set_yticks([])
		grid[8-i].set_xticks([])		





def rho_arrow_plot(dirname,rho_names,vel_names,pixels,color):
	'''probably going to delete this function'''
	fig,grid,im=cube_plot(dirname,rho_names,'rho',pixels)
	add_arrows(grid,dirname,vel_names,int(pixels/40),25,color) #35 for mu=5, 25 for mu=20	

def rho_quiver_plot(dirname,rho_names,vel_names,pixels,spacing,crop):
	'''density imshow with velocity arrows on top'''
	fig,grid,im=cube_plot(dirname,rho_names,'rho',pixels)
	print('Arrows:')
	for i in range(len(rho_names)):
		print(i)
		arrows_version2(grid[i],dirname,vel_names[i],spacing,crop)


def velocity_quiver_plot(dirname,vel_names,pixels,spacing,crop):
	'''velocity magnitude imshow with directional quiver over the top'''
	fig,grid,im=cube_plot(dirname,vel_names,'v',pixels)
	print('Arrows:')
	for i in range(len(vel_names)):
		print(i)
		arrows_version2(grid[i],dirname,vel_names[i],spacing,crop)
		



'''//////////Slicing through grid projections/////////////'''

def spacial_average(name,zoomzone,pixels):
	'''slice through central z axis gives xy plane, average the B_rot/B_pol ratio for each R(x,y)'''
	
	bratio=np.absolute(ratioB_cube(name,zoomzone,pixels))
	mid=int(pixels/2)
	width=int(pixels/40)
	bratio=bratio[:,:,mid-width:mid+width]
	bratio=np.sum(bratio,2).flatten()
	bratio=bratio/(2*width)
	x=np.linspace(-zoomzone,zoomzone,pixels)
	y,x=np.meshgrid(x,x)
	x=x.flatten()
	y=y.flatten()
	
	rs=np.sqrt(y**2+x**2)*code_units.d_cu*1e-17 #graph wants weird units
	mask=np.where(rs<=0.05)
	rs=rs[mask]
	bratio=bratio[mask]
	d=avline(rs,bratio,20)

	return d[1][:-1],d[0]
	
def spacial_multiplot(dirname,snaps,zoomzone,pixels):
	fig,ax=plt.subplots(1)
	fig.set_label('r/ 10^17 cm')
	fig.set_label('B_rot/B_pol')
	for i in range(len(snaps)):
		print(i)
		x,y=spacial_average(dirname+snaps[i],zoomzone,pixels)
		ax.plot(x,y,label=snaps[i][-3:])
	ax.semilogy()
	ax.set_xlim(0,0.05)
	ax.legend()


'''//////////Checking the divergence of the B field/////////////'''

def divB(name):
	a=arepo_utils.aread(name)
	div=np.absolute(a.divb)
	scale=(a.mass/a.rho)**(1./3.)
	mag=np.sqrt(a.bfield[:,0]**2+a.bfield[:,1]**2+a.bfield[:,2]**2)
	divb=div*scale/mag
	return divb,mag





'''//////////Finding best output snapshots times for higher resolution repeats/////////////'''

def track_max_rho(dirname):
	'''plots the maximum denity of snapshots within a directory as a function of time'''
	
	names=np.asarray(glob.glob(dirname+'/snapshot_*'))
	NAMES=[]
	for i in range(len(names)): #extract snap numbers from the snap names 
		N=names[i].rfind('_')
		number=names[i][N+1:]
		NAMES=np.append(NAMES,number)

	args=np.asarray(NAMES.argsort()).astype(int) #sort names in order of snap number 
	rho=np.array([])
	t=np.array([])

	text_trap = io.StringIO() #prevent massive text output from snapshot reads
	sys.stdout = text_trap
	for i in range (len(args)):
		name=names[args[i]]
		a=arepo_utils.aread(name)
		rho=np.append(rho,a.rho.max())
		t=np.append(t,a.time)
	sys.stdout = sys.__stdout__
	#plt.plot(t,rho)
	#plt.plot(t,rho,'x')
	return t,rho

def output_times(dirname,detailed_interval):
	'''takes the maximum desnity time evolution from track_max_rho()
	and fits a spline to it in log space, creates a TIMES.txt file with
	snapshot creation times that undersamples the early (boring) stages of 
	the simulation, and linearly samples the later stages of collapse'''
	t,rho=track_max_rho(dirname)
	spline=interp1d((t),(rho),bounds_error=False,fill_value='extrapolate')
	tnew=np.linspace(t[0],t[-1],1000)
	rhonew=spline(tnew)
	
	plt.figure(),plt.plot(tnew,rhonew,label='low resolution run',color='r')

	grads=np.array([0])
	steepening=np.array([])
	for i in range (len(tnew)-2):
		dy=np.log10(rhonew[i+2])-np.log10(rhonew[i+1])
		dx=np.log10(tnew[i+2])-np.log10(tnew[i+1])
		grads=np.append(grads,dy/dx)
		if dy/np.log10(rhonew[i+1])>=100:
			steepening=np.append(steepening,tnew[i+1])

	mask=np.where(grads<100)
	if len(mask[0])==len(grads):
		print('no density steepening detected')
		time_interesting=tnew.max()*2/3
	else:
		time_interesting=tnew[mask].max()
		#time_interesting=steepening.min()
		print('density steepending detected at t= %.3f'%time_interesting)

	plt.axvline(x=time_interesting,label='begin linear regime')
	plt.legend()

	#tmax=tnew.max()
	mask_sink=np.where(rhonew>2.0e8)
	if len(mask_sink[0])==0:
		print('no sink found, using latest snapshot')
		t_sink=tnew.max()
		tmax=tnew.max()+20*detailed_interval
	else:
		t_sink=tnew[mask_sink].min()
		tmax=t_sink+20*detailed_interval
		print('sink found at t=%.3f'%t_sink)

	#TIMES=np.log(np.linspace(np.exp(t[0]),np.exp(time_interesting),10))
	TIMES=np.log(np.linspace(np.exp(t[0]),np.exp(time_interesting),10))
	TIMES=np.append(TIMES,np.arange(time_interesting,tmax,detailed_interval))

	plt.figure()
	for i in range(len(TIMES)):	
		if i==0:
			plt.axvline(x=TIMES[i]*code_units.t_cu/(60*60*24*365),linewidth=1,label='output times')
		else:
			plt.axvline(x=TIMES[i]*code_units.t_cu/(60*60*24*365),linewidth=1)	
	#plt.scatter(((TIMES)*code_units.t_cu),((spline(TIMES))*code_units.rho_cu),s=0.1,c='r',label='selected outputs')
	plt.plot((tnew*code_units.t_cu/(60*60*24*365)),(rhonew*code_units.rho_cu),label='low resolution run',color='r')
	plt.axhline(2e8*code_units.rho_cu,label='Sink particle creation',color='g')
	plt.yscale('log')
	
	plt.ylabel('log10(rho/gcm^-3)'),plt.xlabel('t/yrs')	
	plt.legend()
	
	f= open(dirname+"TIMES.txt","w+")
	for i in range (len(TIMES)):
		f.write('%s\r\n'%TIMES[i])
	return TIMES,spline(TIMES)





'''//////////tracking the mass of sink particles/////////////'''

def sink_track(dirname):
	names=np.asarray(glob.glob(dirname+'/snapshot_*'))
	NAMES=[]
	for i in range(len(names)): #extract snap numbers from the snap names
		N=names[i].rfind('_')
		number=names[i][N+1:]
		NAMES=np.append(NAMES,number)
	args=np.asarray(NAMES.argsort()).astype(int) #sort names in order of snap number
	M=np.array([])
	t=np.array([])
	N=0

	a=arepo_utils.aread(names[args[-1]]) #check end number of sinks 
	N_sink_end=a.npart[5]
	M=np.zeros((N_sink_end,len(args)))

	text_trap = io.StringIO() #prevent massive text output from snapshot reads
	sys.stdout = text_trap
	for i in range(len(args)):
		a=arepo_utils.aread(names[args[i]])
		Nsink=a.npart[5]
		t=np.append(t,a.time*code_units.t_cu/(60*60*24*365))
		if Nsink>0:
			M[0:Nsink,i]=a.sinkmass
				
			
	sys.stdout = sys.__stdout__
	
	mask=np.where(M[0,:]>0)
	t=t[mask]
	M=M[:,mask]
	plt.figure(),plt.title('first sink at %.0f yrs'%t[0])
	
	for i in range(N_sink_end):
		plt.plot(t,M[i,:][0]),plt.xlabel('time/yrs'),plt.ylabel('number of sinks')
		
	
	return t,M

