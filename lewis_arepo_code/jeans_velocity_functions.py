import numpy as np
import matplotlib.pyplot as plt
import arepo_utils 
import io
import sys
import code_units

def solenoidal_velocity(a):
	'''calculates radial and non radial component of the velocity, also returns total vel magnitude''' 
	#create vector to center 
	mid=np.where(a.rho==a.rho.max())
	xmid=a.x[mid]
	ymid=a.y[mid]
	zmid=a.z[mid]
	r=np.array([a.x-xmid, a.y-ymid, a.z-zmid])
	#radial and non-radial velocities 
	dot_product=(r[0]*a.vx + r[1]*a.vy + r[2]*a.vz)
	cross_product=(   (r[1]*a.vz - r[2]*a.vy) , - (r[0]*a.vz - r[2]*a.vx) , (r[0]*a.vy - r[1]*a.vx)   )
	
	#normalise 
	R=np.sqrt(r[0]**2+r[1]**2+r[2]**2)
	radial_v=abs(dot_product/R)
	nonrad_v=abs(cross_product/R)
	#remove NaNs
	radial_v[np.where(R==0)]=0
	nonrad_v[:,np.where(R==0)]=np.zeros((3,1,1))
	#reduce cross product to 1D
	nonrad_v=np.sqrt(nonrad_v[0]**2 + nonrad_v[1]**2 + nonrad_v[2]**2) 

	v_mag=np.sqrt(a.vx**2+a.vy**2+a.vz**2)
	return radial_v,nonrad_v,v_mag,R


def equal_volumes_average(v,volumes):
	'''The volumes of the cells are unequal, so a reliable mean cannot be taken.
	By using the smallest cell volume as a base unit, larger cells are split based on their 
	volume. e.g. a cell with volume=10*V_smallest will be replaced with 10 cells, each with the 
	original velocity.
	Using the number of cells to replace each cell with, a weighted average is calculated.'''
	
	
	vol_unit=volumes.min()
	num_cells=(np.round(volumes/vol_unit,0)).astype(int)
	weighted_av=sum(v*num_cells)/sum(num_cells)
	return weighted_av,num_cells
	

def snapname(start,i,interval):
	'''creates snapshot id'''
	n='00'+str(start+i*interval)
	if start+i*interval>9:
		n='0'+str(start+i*interval)
	if start+i*interval>99:
		n=str(start+i*interval)
	return n

def cycle(dirname,start,end,interval,zoomzone,name):

	#create txt file to write data
	f = open(name, "x")
	

	N=int((end-start)/interval)
	Enonrad=[]
	Eradial=[]
	Etot=[]
	t=[]
	nonradial_av=[]
	radial_av=[]
	v_av=[]
	for i in range(N):
		
		text_trap = io.StringIO() #prevent massive text output from snapshot reads
		sys.stdout = text_trap
		#read snap
		n=snapname(start,i,interval)
		a=arepo_utils.aread(dirname+'snapshot_'+n)
		#radial and non-radial velocities 
		radial_v,nonrad_v,v,Rs=solenoidal_velocity(a)
		#equalise volumes and take weighted average
#		radial_v=equal_volumes_average(radial_v, a.mass/a.rho)
#		nonrad_v=equal_volumes_average(nonrad_v, a.mass/a.rho)
#		v=equal_volumes_average(v, a.mass/a.rho)

		mask=np.where(Rs<zoomzone)
		Enonrad_=sum(0.5*a.mass[mask]*nonrad_v[mask]**2)
		Eradial_=sum(0.5*a.mass[mask]*radial_v[mask]**2)
		Etot_=sum(0.5*a.mass[mask]*v[mask]**2)

		nonradial_av_,w=equal_volumes_average(nonrad_v[mask],(a.mass/a.rho)[mask])
		radial_av_,w=equal_volumes_average(radial_v[mask],(a.mass/a.rho)[mask])
		v_av_,w=equal_volumes_average(v[mask],(a.mass/a.rho)[mask])

		#add to time laps arrays
		Enonrad.append(Enonrad_)
		Eradial.append(Eradial_)
		Etot.append(Etot_)
		t.append(a.time)
		nonradial_av.append(nonradial_av_)
		radial_av.append(radial_av_)
		v_av.append(v_av_)


		sys.stdout = sys.__stdout__
		f.write(str(Enonrad_) + ' ' + str(Eradial_) + ' ' + str(Etot_) + ' ' + str(nonradial_av_) + ' ' + str(radial_av_) + ' ' + str(v_av_) + ' ' + str(a.time) + '\n')
		print(n + ' :written')
	return Enonrad_,Eradial_,Etot_,nonradial_av_,radial_av_,v_av_,a.time
		
def txtread(txtfile):
	Enonrad=[]
	Eradial=[]	
	Etot=[]
	vnr=[]
	vr=[]
	v=[]
	t=[]
	with open(txtfile) as f:
		for line in f.readlines():
			Enonrad.append(line.split()[0])
			Eradial.append(line.split()[1])
			Etot.append(line.split()[2])
			vnr.append(line.split()[3])
			vnr.append(line.split()[4])
			vnr.append(line.split()[5])
			t.append(line.split()[6])
	return Enonrad,Eradial,Etot ,vnr,vr,v,t



def weighted_hist(snap,zoomzone,bins):
	'''create PDF of radial, non-radial and total velcoity for a single snapshot'''
	#read and zoom in
	a=arepo_utils.aread(snap)
	radial_v,nonrad_v,v,Rs=solenoidal_velocity(a)
	mask=np.where(Rs<zoomzone)

	#find how many volume units there are per cell for weighted histogram
	average_nonrad,weight_nrad=equal_volumes_average(nonrad_v[mask],(a.mass/a.rho)[mask])
	average_rad,weight_rad=equal_volumes_average(radial_v[mask],(a.mass/a.rho)[mask])
	average_v,weight_v=equal_volumes_average(v[mask],(a.mass/a.rho)[mask])

	#create PDFs
	Av,Bv=np.histogram(v[mask],weights=weight_v,bins=bins,density=True)
	Anv,Bnv=np.histogram(nonrad_v[mask],weights=weight_v,bins=bins,density=True)
	Arv,Brv=np.histogram(radial_v[mask],weights=weight_v,bins=bins,density=True)

	return Av,Bv,Anv,Bnv,Arv,Brv,average_nonrad,average_rad,average_v


	


def PDF_plotter(snaps,zoomzone,bins,labels):
	'''
	snaps='16jeans','32jeans','64jeans','128jeans'
	zoomzone=0.01
	bins=100
	labels=np.array(['16 cells','32 cells','64 cells','128 cells'])
	'''		
	colors=np.array(['b','r','chartreuse','k'])
	fig,axs=plt.subplots(3,sharex=True)
	for i in range(len(snaps)):
		Av,Bv,Anv,Bnv,Arv,Brv,average_nonrad,average_rad,average_v=weighted_hist(snaps[i],zoomzone,bins)
		axs[0].plot(Bv[:-1]*code_units.v_cu/1e5,Av,label=labels[i],color=colors[i],ls='steps')
		axs[1].plot(Bnv[:-1]*code_units.v_cu/1e5,Anv,color=colors[i],ls='steps')
		axs[2].plot(Brv[:-1]*code_units.v_cu/1e5,Arv,color=colors[i],ls='steps') 
		
		
		#axs[0].axvline(average_v*code_units.v_cu/1e5,0,0.15,color=colors[i],linestyle='--')
		#axs[1].axvline(average_nonrad*code_units.v_cu/1e5,0,0.15,color=colors[i],linestyle='--')
		#axs[2].axvline(average_rad*code_units.v_cu/1e5,0,0.15,color=colors[i],linestyle='--')
	axs[0].yaxis.set_visible(False)
	axs[1].yaxis.set_visible(False)
	axs[2].yaxis.set_visible(False)
	axs[1].text(0.8,0.6,r'$\frac{v \times r}{|\ r\ |}$',transform=axs[1].transAxes,fontsize=20)
	axs[2].text(0.8,0.6,r'$\frac{v \cdot r}{|\ r\ |}$',transform=axs[2].transAxes,fontsize=20)
	axs[2].set_xlabel(r'$v \ [kms^{-1}$]',fontsize=15)	
	axs[0].legend(fontsize=11)
	axs[0].tick_params(axis="x", labelsize=15)
	axs[1].tick_params(axis="x", labelsize=15)
	axs[2].tick_params(axis="x", labelsize=15)
	plt.subplots_adjust(wspace=0, hspace=0)
	plt.xlim(0,Brv.max()*code_units.v_cu/1e5)

def average_fromPDF(snaps,zoomzone):
	res=np.array([16,32,64,128])
	rv=[]
	nv=[]
	V=[]
	for i in range(len(snaps)):
		a=arepo_utils.aread(snaps[i])
		radial_v,nonrad_v,v,Rs=solenoidal_velocity(a)
		mask=np.where(Rs<zoomzone)
		avnrad,weight=equal_volumes_average(nonrad_v[mask],(a.mass/a.rho)[mask])
		avrad,weight=equal_volumes_average(radial_v[mask],(a.mass/a.rho)[mask])
		avv,weight=equal_volumes_average(v[mask],(a.mass/a.rho)[mask])
		
		rv.append(avrad)
		nv.append(avnrad)
		V.append(avv)

	fig,ax=plt.subplots(1)
	ax.plot(res,np.asarray(V)*code_units.v_cu/1e5,color='k',label='total velocity')
	ax.plot(res,np.asarray(rv)*code_units.v_cu/1e5,color='r',label='radial component')
	ax.plot(res,np.asarray(nv)*code_units.v_cu/1e5,color='b',label='non-radial component')
	ax.legend(frameon=False)
		
