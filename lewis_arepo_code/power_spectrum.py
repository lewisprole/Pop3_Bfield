from projection_functions import * 
from code_units import * 
from scipy.stats import binned_statistic



fig,ax=plt.subplots(ncols=2,sharey=True)
plt.subplots_adjust(wspace=0)
k,P=create_spectrum_zeropad('/scratch/c.c1521474/resolution_test/MHD/1e10MHD/magnetic_grid_001',10)
ax[0].loglog(0.1/k *code_units.d_cu/ap.au.cgs.value,P,'.',label='-250yr')
k,P=create_spectrum_zeropad('/scratch/c.c1521474/resolution_test/MHD/1e10MHD/magnetic_grid_014',10)
ax[0].loglog(0.1/k *code_units.d_cu/ap.au.cgs.value,P,'.',label='0yr')
k,P=create_spectrum_zeropad('/scratch/c.c1521474/resolution_test/MHD/1e10MHD/magnetic_grid_071',10)
ax[0].loglog(0.1/k *code_units.d_cu/ap.au.cgs.value,P,'.',label='250yr')
k,P=create_spectrum_zeropad('/scratch/c.c1521474/resolution_test/MHD/1e10MHD/magnetic_grid_128',10)
ax[0].loglog(0.1/k *code_units.d_cu/ap.au.cgs.value,P,'.',label='500yr')
#k,P=create_spectrum('/scratch/c.c1521474/resolution_test/MHD/1e10MHD/magnetic_grid_242',1,'no')
#ax[0].loglog(0.1/k *code_units.d_cu/ap.au.cgs.value,P,label='1000yr')
ax[0].text(0.7,0.95,r'Kazantsev spectrum',ha='center', va='center', transform=ax[0].transAxes,fontsize=9)
ax[0].legend(frameon=False,loc='lower left',fontsize=8)
k,P=create_spectrum_zeropad('/scratch/c.c1521474/resolution_test/MHD/1e10_uniform/magnetic_grid_001',10)
ax[1].loglog(0.1/k *code_units.d_cu/ap.au.cgs.value,P,'.',label='-250yr')
k,P=create_spectrum_zeropad('/scratch/c.c1521474/resolution_test/MHD/1e10_uniform/magnetic_grid_013',10)
ax[1].loglog(0.1/k *code_units.d_cu/ap.au.cgs.value,P,'.',label='0yr')
k,P=create_spectrum_zeropad('/scratch/c.c1521474/resolution_test/MHD/1e10_uniform/magnetic_grid_071',10)
ax[1].loglog(0.1/k *code_units.d_cu/ap.au.cgs.value,P,'.',label='250yr')
k,P=create_spectrum_zeropad('/scratch/c.c1521474/resolution_test/MHD/1e10_uniform/magnetic_grid_128',10)
ax[1].loglog(0.1/k *code_units.d_cu/ap.au.cgs.value,P,'.',label='500yr')
ax[1].text(0.8,0.95,r'Uniform field',ha='center', va='center', transform=ax[1].transAxes,fontsize=9)
#ax[1].legend(frameon=False,loc='lower left',fontsize=8)
ax[0].tick_params(axis="x", labelsize=10,direction="in",which='both')
ax[0].tick_params(axis="y", labelsize=10,direction="in",which='both')
ax[1].tick_params(axis="y", labelsize=10,direction="in",which='both')
ax[1].tick_params(axis="x", labelsize=10,direction="in",which='both')
ax[0].set_xlim(0.11/1 * d_cu/ap.au.cgs.value, 0.1/500 * d_cu/ap.au.cgs.value)
ax[1].set_xlim(0.11/1 * d_cu/ap.au.cgs.value, 0.1/500 * d_cu/ap.au.cgs.value)
ax[1].set_xlabel('$\lambda$ [AU]',fontsize=10)
ax[0].set_xlabel('$\lambda$ [AU]',fontsize=10)
ax[0].set_ylabel(r'P$_{\rm B}$        ',fontsize=10,rotation=0)




fig,ax=plt.subplots(nrows=2,sharex=True)
plt.subplots_adjust(hspace=0)
rhos=r'$\rho_{sink}$=10$^{-10}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-9}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-8}$gcm$^{-3}$'
k,P=create_spectrum_zeropad('/scratch/c.c1521474/resolution_test/MHD2/1e8MHD/magnetic_grid_009',150)
ax[0].loglog(0.002/(100/(400/k))*d_cu/ap.au.cgs.value,P,'x')
k,P=create_spectrum_zeropad('/scratch/c.c1521474/resolution_test/MHD2/1e9MHD/magnetic_grid_051',150)
ax[0].loglog(0.002/(100/(400/k))*d_cu/ap.au.cgs.value,P,'x')
k,P=create_spectrum_zeropad('/scratch/c.c1521474/resolution_test/MHD2/1e10MHD/magnetic_grid_044',150)
ax[0].loglog(0.002/(100/(400/k))*d_cu/ap.au.cgs.value,P,'x')
ax[0].set_xlim(620,2.673834848907378)
ax[0].set_ylabel(r'P$_{\rm B}                $',fontsize=10,rotation=0)
ax[1].set_ylabel(r'P$_{\rm B}                $',fontsize=10,rotation=0)
ax[1].set_xlabel('$\lambda$ [AU]',fontsize=10)
ax[0].loglog(0.002/(100/(400/k))*d_cu/ap.au.cgs.value,((100/(400/k)))**(3/2)*4e5,'--',c='k')
k,P=create_spectrum_zeropad('/scratch/c.c1521474/resolution_test/MHD2/1e8_uniform/magnetic_grid_009',100)
ax[1].loglog(0.002/(100/(400/k))*d_cu/ap.au.cgs.value,P,'x',label=rhos[0])
k,P=create_spectrum_zeropad('/scratch/c.c1521474/resolution_test/MHD2/1e9_uniform/magnetic_grid_051',100)
ax[1].loglog(0.002/(100/(400/k))*d_cu/ap.au.cgs.value,P,'x',label=rhos[1])
k,P=create_spectrum_zeropad('/scratch/c.c1521474/resolution_test/MHD2/1e10_uniform/magnetic_grid_043',100)
ax[1].loglog(0.002/(100/(400/k))*d_cu/ap.au.cgs.value,P,'x',label=rhos[2])
ax[1].legend(frameon=False,loc='lower left',fontsize=8)
ax[0].text(0.1,0.5,r'$\propto k^{3/2}$',ha='center', va='center', transform=ax[0].transAxes,fontsize=10)
ax[0].text(0.85,0.95,r'Kazantsev spectrum',ha='center', va='center', transform=ax[0].transAxes,fontsize=10)
ax[1].text(0.9,-0.05,r'Uniform field',ha='center', va='center', transform=ax[0].transAxes,fontsize=10)
ax[0].tick_params(axis="x", labelsize=10,direction="in",which='both')
ax[1].tick_params(axis="y", labelsize=10,direction="in",which='both')


fig,ax=plt.subplots(nrows=2,sharex=True)
plt.subplots_adjust(hspace=0)
files='/scratch/c.c1521474/resolution_test/MHD2/1e8MHD/snapshot_009','/scratch/c.c1521474/resolution_test/MHD2/1e9MHD/snapshot_051','/scratch/c.c1521474/resolution_test/MHD2/1e10MHD/snapshot_044'
files2='/scratch/c.c1521474/resolution_test/MHD2/1e8_uniform/snapshot_009','/scratch/c.c1521474/resolution_test/MHD2/1e9_uniform/snapshot_051','/scratch/c.c1521474/resolution_test/MHD2/1e10_uniform/snapshot_044'
rhos=r'$\rho_{sink}$=10$^{-10}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-9}$gcm$^{-3}$',r'$\rho_{sink}$=10$^{-8}$gcm$^{-3}$'
for i in range(len(files)):
	a=aread(files[i])
	midx,midy,midz=a.sinkx[0],a.sinky[0],a.sinkz[0]
	r=np.sqrt((a.x-midx)**2+(a.y-midy)**2+(a.z-midz)**2)
	b,rs,y=binned_statistic(r,a.bmag,bins=10**np.linspace(np.log10(r.min()),np.log10(r.max()),25))
	ax[0].loglog(rs[:-1]*code_units.d_cu /ap.au.cgs.value ,b*code_units.B_cu)
	a=aread(files2[i])
	midx,midy,midz=a.sinkx[0],a.sinky[0],a.sinkz[0]
	r=np.sqrt((a.x-midx)**2+(a.y-midy)**2+(a.z-midz)**2)
	b,rs,y=binned_statistic(r,a.bmag,bins=10**np.linspace(np.log10(r.min()),np.log10(r.max()),25))
	ax[1].loglog(rs[:-1]*code_units.d_cu /ap.au.cgs.value ,b*code_units.B_cu)



