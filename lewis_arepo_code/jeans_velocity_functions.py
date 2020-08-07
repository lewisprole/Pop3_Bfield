import numpy as np
import matplotlib.pyplot as plt
import arepo_utils 
import io
import sys

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
	return radial_v,nonrad_v,v_mag


def equal_volumes_average(v,volumes):
	'''The volumes of the cells are unequal, so a reliable mean cannot be taken.
	By using the smallest cell volume as a base unit, larger cells are split based on their 
	volume. e.g. a cell with volume=10*V_smallest will be replaced with 10 cells, each with the 
	original velocity.
	Using the number of cells to replace each cell with, a weighted average is calculated.'''
	
	
	vol_unit=volumes.min()
	num_cells=(np.round(volumes/vol_unit,0)).astype(int)
	weighted_av=sum(v*num_cells)/sum(num_cells)
	return weighted_av
	

def snapname(start,i,interval):
	'''creates snapshot id'''
	n='00'+str(start+i*interval)
	if start+i*interval>9:
		n='0'+str(start+i*interval)
	if start+i*interval>99:
		n=str(start+i*interval)
	return n

def cycle(dirname,start,end,interval,name):

	#create txt file to write data
	f = open(name, "x")
	

	N=int((end-start)/interval)
	Enonrad=[]
	Eradial=[]
	Etot=[]
	t=[]
	for i in range(N):
		
		text_trap = io.StringIO() #prevent massive text output from snapshot reads
		sys.stdout = text_trap
		#read snap
		n=snapname(start,i,interval)
		a=arepo_utils.aread(dirname+'snapshot_'+n)
		#radial and non-radial velocities 
		radial_v,nonrad_v,v=solenoidal_velocity(a)
		#equalise volumes and take weighted average
#		radial_v=equal_volumes_average(radial_v, a.mass/a.rho)
#		nonrad_v=equal_volumes_average(nonrad_v, a.mass/a.rho)
#		v=equal_volumes_average(v, a.mass/a.rho)


		Enonrad_=sum(0.5*a.mass*nonrad_v**2)
		Eradial_=sum(0.5*a.mass*radial_v**2)
		Etot_=sum(0.5*a.mass*v**2)

		#add to time laps arrays
		Enonrad.append(Enonrad_)
		Eradial.append(Eradial_)
		Etot.append(Etot_)
		t.append(a.time)

		sys.stdout = sys.__stdout__
		f.write(str(Enonrad_) + ' ' + str(Eradial_) + ' ' + str(Etot_) + ' ' + str(a.time) + '\n')
		print(n + ' :written')
	return Enonrad_,Eradial_,Etot_,a.time 
		
def txtread(txtfile):
	Enonrad=[]
	Eradial=[]	
	Etot=[]
	t=[]
	with open(txtfile) as f:
		for line in f.readlines():
			Enonrad.append(line.split()[0])
			Eradial.append(line.split()[1])
			Etot.append(line.split()[2])
			t.append(line.split()[3])
	return Enonrad,Eradial,Etot ,t
