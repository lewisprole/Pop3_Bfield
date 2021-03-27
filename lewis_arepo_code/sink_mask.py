import numpy as np 




def sinkmask(a,zoomzone):
	midx=a.sinkx[np.where(a.sinkmass==a.sinkmass.max())]
	midy=a.sinky[np.where(a.sinkmass==a.sinkmass.max())]
	midz=a.sinkz[np.where(a.sinkmass==a.sinkmass.max())]
	r=np.sqrt((a.x-midx)**2+(a.y-midy)**2+(a.z-midz)**2)
	mask=np.where(r<zoomzone)
	return mask
