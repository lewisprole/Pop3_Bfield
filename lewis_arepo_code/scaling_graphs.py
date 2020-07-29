import numpy as np
import matplotlib.pyplot as plt 
import os.path, time

def speedup(dir_main,dirs_specific,start_no,end_no):
	T=np.array([])
	for i in range(len(dirs_specific)):
		t0=time.ctime(os.path.getctime(dir_main+dirs_specific[i]+'snapshot_'+start_no))
		t1=time.ctime(os.path.getctime(dir_main+dirs_specific[i]+'snapshot_'+end_no))

		t0=t0[-13:-5]
		t1=t1[-13:-5]

		hdiff=int(t1[:2])-int(t0[:2])
		mdiff=int(t1[3:5])-int(t0[3:5])
		sdiff=int(t1[6:8])-int(t0[6:8])
		diff=hdiff*60*60 + mdiff*60 + sdiff
		T=np.append(T,diff)
	return T 






