import arepo_utils 
from matplotlib.patches import Circle
from matplotlib.lines import Line2D



def plot(file)
	a=aread(file)
	mid=np.where(a.rho==a.rho.max())
	r=np.sqrt((a.x-a.x[mid])**2+(a.y-a.y[mid])**2+(a.z-a.z[mid])**2)
	heat=np.array([8,9,10,11,12,20,21,22,23,24,25])
	cool=np.array([0,1,2,3,4,6,7,13,14,15,16,17,18,19,26])
	tff=np.sqrt(3*np.pi/(32*ap.G.cgs.value*a.rho*rho_cu))
	x,y,z=np.histogram2d(np.log10(abs(a.u/np.sum(a.cooling[:,heat],1))*t_cu/(60*60*24*365)),np.log10(a.rho*rho_cu),bins=(800,800))
        line1=plt.imshow(x/x,cmap='autumn',aspect='auto',label=r'$t_h$',extent=[z[0],z[-1],y[-1],y[0]])
	x,y,z=np.histogram2d(np.log10(abs(a.u/np.sum(a.cooling[:,cool],1))*t_cu/(60*60*24*365)),np.log10(a.rho*rho_cu),bins=(800,800))
	line2=plt.imshow(x/x,cmap='winter',aspect='auto',label=r'$t_c$',extent=[z[0],z[-1],y[-1],y[0]])
	plt.ylim(y[0],y[-1])
	rho=10**np.linspace(z[0],z[-1],100)
	tff=np.sqrt(3*np.pi/(32*ap.G.cgs.value*rho))/(60*60*24*365)
	line3=plt.plot(np.log10(rho),np.log10(tff),'k',label=r'$t_{ff}$')
	line1=Line2D([0], [0], color='r', lw=2)
	line2=Line2D([0], [0], color='b', lw=2)
	line3=Line2D([0], [0], color='k', lw=2)
	plt.legend([line1,line2,line3],(r'$t_h$',r'$t_c$',r'$t_{ff}$'),fontsize=12,frameon=False,markerscale=10)
	plt.subplots_adjust(left = 0.2,bottom = 0.17,right=0.9)
        plt.tick_params(axis="x", labelsize=15,direction="in")
        plt.tick_params(axis="y", labelsize=15,direction="in")
        plt.xlabel(r'log$_{10}(\rho$ [gcm$^{-3}$])',fontsize=20)
        plt.ylabel(r'log$_{10}(t$ [yrs])       ',fontsize=20,rotation=90)
	plt.ylim(np.log10(tff).min()-1,y[-1])




