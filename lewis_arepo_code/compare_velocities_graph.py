import numpy as np
import matplotlib.pyplot as plt 
import arepo_utils 
import code_units 
import projection_functions 
import read_sink_info
from mpl_toolkits.axes_grid1 import make_axes_locatable

#set up figure space 

fig = plt.figure()
ax1a=plt.subplot2grid((3, 8), (0, 0), colspan=1,rowspan=1)
ax1b=plt.subplot2grid((3, 8), (0, 1), colspan=1,rowspan=1)
ax1c=plt.subplot2grid((3, 8), (1, 0), colspan=2,rowspan=1)
ax1d=plt.subplot2grid((3, 8), (2, 0), colspan=2,rowspan=1)#,sharex=ax1c)

ax1a.set_yticks([])
ax1b.set_yticks([])
ax1a.set_xticks([])
ax1b.set_xticks([])

ax2a=plt.subplot2grid((3, 8), (0, 3), colspan=1,rowspan=1)
ax2b=plt.subplot2grid((3, 8), (0, 4), colspan=1,rowspan=1)
ax2c=plt.subplot2grid((3, 8), (1, 3), colspan=2,rowspan=1)
ax2d=plt.subplot2grid((3, 8), (2, 3), colspan=2,rowspan=1)#,sharex=ax2c)

ax2a.set_yticks([])
ax2b.set_yticks([])
ax2a.set_xticks([])
ax2b.set_xticks([])

ax3a=plt.subplot2grid((3, 8), (0, 6), colspan=1,rowspan=1)
ax3b=plt.subplot2grid((3, 8), (0, 7), colspan=1,rowspan=1)
ax3c=plt.subplot2grid((3, 8), (1, 6), colspan=2,rowspan=1)
ax3d=plt.subplot2grid((3, 8), (2, 6), colspan=2,rowspan=1)#,sharex=ax3c)

ax3a.set_yticks([])
ax3b.set_yticks([])
ax3a.set_xticks([])
ax3b.set_xticks([])

plt.subplots_adjust(wspace=0, hspace=0,bottom=0.15,top=0.95)

#labels

ax1c.set_ylabel(r'N$_{\rm sink}$',fontsize=10)
ax1d.set_ylabel(r'M$_{\rm tot}$ [M$_\odot$]',fontsize=10)
ax1d.set_xlabel('t [yr]',fontsize=10)

ax2c.set_ylabel(r'N$_{\rm sink}$',fontsize=10)
ax2d.set_ylabel(r'M$_{\rm tot}$ [M$_\odot$]',fontsize=10)
ax2d.set_xlabel('t [yr]',fontsize=10)

ax3c.set_ylabel(r'N$_{\rm sink}$',fontsize=10)
ax3d.set_ylabel(r'M$_{\rm tot}$ [M$_\odot$]',fontsize=10)
ax3d.set_xlabel('t [yr]',fontsize=10)

#plot projections 

dirname='/scratch/c.c1521474/resolution_test/'
#A12=dirname+'/merge/1e12_projections/density_grid_495'
A12=dirname+'seed1_alpha0.25/density_grid_495'
A12_=dirname+'seed1_alpha0.25/density_grid_072'
B12=dirname+'/seed4/1e12/density_grid_488'
C12=dirname+'/seed5/1e12/density_grid_445'

rho=projection_functions.read_cube(A12)
rho=rho*code_units.rho_cu
im=np.log10(np.sum(rho,1))
V1,V2=im.min(),im.max()
im=ax1a.imshow(np.log10(np.sum(rho,1)),vmin=V1,vmax=V2)
x,y,z,M=projection_functions.read_sink(A12[:-16]+'/',A12[-3:])
ax1a.scatter(z,x,s=0.5,c='magenta')
ax1a.set_xlim(0,500)
ax1a.set_ylim(0,500)


rho=projection_functions.read_cube(A12_)
rho=rho*code_units.rho_cu
im=ax1b.imshow(np.log10(np.sum(rho,0)),vmin=V1,vmax=V2)
x,y,z,M=projection_functions.read_sink(A12_[:-16]+'/',A12_[-3:])
ax1b.scatter(z,y,s=0.5,c='magenta')
ax1b.set_xlim(0,500)
ax1b.set_ylim(0,500)
ax1b.set_aspect('auto')
divider = make_axes_locatable(ax1b)
cax = divider.append_axes("right", size="5%", pad=0)
cbar=fig.colorbar(im, cax=cax,ticks=[-11, -10, -9,-8])
cbar.ax.tick_params(labelsize=9)
cbar.set_label(r'log($\rho$ [gcm$^{-3}$])', rotation=270,fontsize=9i,labelpad=15
ax1a.set_aspect('auto')

rho=projection_functions.read_cube(B12)
rho=rho*code_units.rho_cu
im=np.log10(np.sum(rho,2))
V1,V2=im.min(),im.max()
ax2a.imshow(np.log10(np.sum(rho,2)),vmin=V1,vmax=V2)
x,y,z,M=projection_functions.read_sink(B12[:-16]+'/',B12[-3:])
ax2a.scatter(y,x,s=0.5,c='magenta')
ax2a.set_xlim(0,500)
ax2a.set_ylim(0,500)

rho=projection_functions.read_cube(C12)
rho=rho*code_units.rho_cu
im=np.log10(np.sum(rho,2))
V1,V2=im.min(),im.max()
ax3a.imshow(np.log10(np.sum(rho,2)),vmin=V1,vmax=V2)
x,y,z,M=projection_functions.read_sink(C12[:-16]+'/',C12[-3:])
ax3a.scatter(y,x,s=0.5,c='magenta')
ax3a.set_xlim(0,500)
ax3a.set_ylim(0,500)


#plot evolutions 
extensions1= 'merge/1e12//sink_particle_info/','seed4/1e12//sink_particle_info/','seed5/1e12//sink_particle_info/'
extensions2 = '/seed1_alpha0.25//sink_particle_info/','/seed4_alpha0.25//sink_particle_info/','/seed5_alpha0.25//sink_particle_info/'

t,N,M=read_sink_info.Nsinks(dirname+extensions1[0])
ax1c.plot((t-t[0])*code_units.t_cu/(60*60*24*365),N,'b')
ax1d.plot((t-t[0])*code_units.t_cu/(60*60*24*365),M,'b')
t,N,M=read_sink_info.Nsinks(dirname+extensions2[0])
ax1c.plot((t-t[0])*code_units.t_cu/(60*60*24*365),N,'b',linestyle='--')
ax1d.plot((t-t[0])*code_units.t_cu/(60*60*24*365),M,'b',linestyle='--')
ax1d.set_xlim(0,400)
ax1c.set_xlim(0,400)
ax1c.set_xticks([])

t,N,M=read_sink_info.Nsinks(dirname+extensions1[1])
ax2c.plot((t-t[0])*code_units.t_cu/(60*60*24*365),N,'r')
ax2d.plot((t-t[0])*code_units.t_cu/(60*60*24*365),M,'r')
t,N,M=read_sink_info.Nsinks(dirname+extensions2[1])
ax2c.plot((t-t[0])*code_units.t_cu/(60*60*24*365),N,'r',linestyle='--')
ax2d.plot((t-t[0])*code_units.t_cu/(60*60*24*365),M,'r',linestyle='--')
ax2d.set_xlim(0,400)
ax2c.set_xlim(0,400)
ax2c.set_xticks([])

t,N,M=read_sink_info.Nsinks(dirname+extensions1[2])
ax3c.plot((t-t[0])*code_units.t_cu/(60*60*24*365),N,'g')
ax3d.plot((t-t[0])*code_units.t_cu/(60*60*24*365),M,'g')
t,N,M=read_sink_info.Nsinks(dirname+extensions2[2])
ax3c.plot((t-t[0])*code_units.t_cu/(60*60*24*365),N,'g',linestyle='--')
ax3d.plot((t-t[0])*code_units.t_cu/(60*60*24*365),M,'g',linestyle='--')
ax3d.set_xlim(0,400)
ax3c.set_xlim(0,400)
ax3c.set_xticks([])
