
import numpy as np
import matplotlib.pyplot as plt 
import arepo_utils 


def snapname(start,i,interval):
	'''creates snapshot id'''
	n='00'+str(start+i*interval)
	if start+i*interval>9:
		n='0'+str(start+i*interval)
	if start+i*interval>99:
		n=str(start+i*interval)
	return n


def cycle(dirname,start,end,interval):
	Mtot=[]
	N=[]
	t=[]
	num=int((end-start)/interval)
	for i in range(num):
		text_trap = io.StringIO() #prevent massive text output from snapshot reads
		sys.stdout = text_trap
		n=snapname(start,i,interval)
		a=arepo_utils.aread(dirname+'snapshot_'+n)
		if a.npart[-1]>0:
			Mtot.append(sum(a.sinkmass))
			N.append(a.npart[-1])
			t.append(a.time)
		sys.stdout = sys.__stdout__
	return Mtot,N,t		



