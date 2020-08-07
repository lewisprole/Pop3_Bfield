import numpy as np
import arepo_utils 


def snapname(start,i):
	n='00'+str(start+i)
	if start+i>9:	
		n='0'+str(start+i)
	if start+i>99:
		n=str(start+i)
	return n

def cycle(dirname,start,end):
	N=int(end-start)
	for i in range(N):
		#read snap
		n=snapname(start,i)
		a=arepo_utils.aread(dirname+'snapshot_'+n)
		
		
