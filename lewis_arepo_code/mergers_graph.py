import numpy as np 
import projection_functions

def graph(files_merge,snapshot_merge,boxsize,files_nomerge)
	fig,ax=plt.subplots(ncols=2,nrows=2)
	plt.subplots_adjust(hspace=0,wspace=0)
	for i in range(len(files_merge)):
		rho=projection_functions.read_cube(files_merge[i])
		flat=np.sum(rho,2)/len(rho[:,:,0])
		ax[0,i].imshow(np.log10(flat),cmap='magma')
		if i>0:
			a=arepo_utils.aread(snapshot_merge[i])
			x=a.sinkx/boxsize*500
			y=a.sinky/boxsize*500
			ax[0,i].scatter(y,x,s=1,c='g')
		ax[0,i].set_yticks([])
		ax[0,i].set_xticks([])
		rho=resolution_test_functions.read_cube(files_nomerge[i])
		flat=np.sum(rho,2)/len(rho[:,:,0])
		ax[1,i].imshow(np.log10(flat),,cmap='magma')
		ax[1,i].set_yticks([])
		ax[1,i].set_xticks([])
		
	
