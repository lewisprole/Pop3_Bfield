import numpy as np
import matplotlib.pyplot as plt
import projection_functions 
import arepo_utils 
import code_units 

def prep(name):
	data=projection_functions.read_cube(name)
	N=len(data[:,:,0])
	im=np.sum(data,2)/N
	return im

def snapname(no):
	'''creates snapshot id'''
	n='00'+str(no)
	if no>9:
		n='0'+str(no)
	if no>99:
		n=str(no)
	return n

def sink_pos(snap,pixels,imx1,imx2,imy1,imz1):
	a=arepo_utils.aread(snap)
	if a.npart[-1]>0:
		X,Y,Z=a.sinkx-imx1,a.sinky-imy1,a.sinkz-imz1
		pixel=(imx2-imx1)/pixels
		x,y,z=(X/pixel).astype(int),(Y/pixel).astype(int),(Z/pixel).astype(int)
	else:
		x,y,z=np.array([0]),np.array([0]),np.array([0])
	return x,y,z

def plot_save(dirname,imnames,nos,imx1,imx2,imy1,imz1,outdir):
	for i in range(len(nos)):
		fig,ax=plt.subplots(1)
		n=snapname(nos[-i-1])
		im=prep(dirname+imnames+n)
		im=im*code_units.rho_cu
		x,y,z=sink_pos(dirname+'snapshot_'+n,1000,imx1,imx2,imy1,imz1)
		if i==0:
			vmin=(np.log10(im).min())
			vmax=(np.log10(im).max())
		image=ax.imshow(np.log10(im),vmin=vmin,vmax=vmax,cmap='afmhot')
		cbar=plt.colorbar(image,pad=0)
		plt.text(1200,600,r'log$_{10}$( $\rho$ [gcm$^{-3}] )$', rotation=270)
		ax.scatter(y,x,c='royalblue',s=1)
		ax.set_ylim(0,1000)
		ax.set_xlim(0,1000)
		ax.set_aspect('equal')
		ax.set_yticks([])
		ax.set_xticks([])
		plt.savefig(outdir+n+'.png')
		plt.close('all')
	


