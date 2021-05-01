import numpy as np 
import projection_functions
import arepo_utils
import matplotlib.pyplot as plt  

def graph(files_merge1,snapshot_merge1,files_nomerge1,snapshot_nomerge1,files_merge2,snapshot_merge2,files_nomerge2,snapshot_nomerge2,boxsize):
	fig = plt.figure()
	ax1=plt.subplot2grid((13, 9), (0, 0), colspan=3,rowspan=3)
	ax2=plt.subplot2grid((13, 9), (0, 3), colspan=3,rowspan=3)
	ax3=plt.subplot2grid((13, 9), (0, 6), colspan=3,rowspan=3)
	ax4=plt.subplot2grid((13, 9), (3, 0), colspan=3,rowspan=3)
	ax5=plt.subplot2grid((13, 9), (3, 3), colspan=3,rowspan=3)
	ax6=plt.subplot2grid((13, 9), (3, 6), colspan=3,rowspan=3)

	ax7=plt.subplot2grid((13, 9), (7, 0), colspan=3,rowspan=3)
	ax8=plt.subplot2grid((13, 9), (7, 3), colspan=3,rowspan=3)
	ax9=plt.subplot2grid((13, 9), (7, 6), colspan=3,rowspan=3)
	ax10=plt.subplot2grid((13, 9), (10, 0), colspan=3,rowspan=3)
	ax11=plt.subplot2grid((13, 9), (10, 3), colspan=3,rowspan=3)
	ax12=plt.subplot2grid((13, 9), (10, 6), colspan=3,rowspan=3)

#	for i in range(3):
#		ax[2,i].spines["left"].set_visible(False)
#		ax[2,i].spines["right"].set_visible(False)
#		ax[2,i].set_yticks([])
	for i in range(len(files_merge1)):
		for j in range(2):
			if j==0:
				files_merge=files_merge1
				ax_merge=np.array([ax1,ax2,ax3])[-i-1]
				snapshot_merge=snapshot_merge1
				files_nomerge=files_nomerge1
				ax_nomerge=np.array([ax4,ax5,ax6])[-i-1]
				snapshot_nomerge=snapshot_nomerge1
			else:
				files_merge=files_merge2
				ax_merge=np.array([ax7,ax8,ax9])[-i-1]
				snapshot_merge=snapshot_merge2
				files_nomerge=files_nomerge2
				ax_nomerge=np.array([ax10,ax11,ax12])[-i-1]
				snapshot_nomerge=snapshot_nomerge2

			rho=projection_functions.read_cube(files_merge[-i-1])
			flat=np.sum(rho,2)/len(rho[:,:,0])
			if i==0:
				vmin=np.log10(flat).min()
				vmax=np.log10(flat).max()

			ax_merge.imshow(np.log10(flat),cmap='magma',vmin=vmin,vmax=vmax)
			if i<2:
				a=arepo_utils.aread(snapshot_merge[-i-1])
				x=a.sinkx/boxsize*500
				y=a.sinky/boxsize*500
				ax_merge.scatter(y,x,s=1,c='g')
			
				a=arepo_utils.aread(snapshot_nomerge[-i-1])
				x=a.sinkx/boxsize*500
				y=a.sinky/boxsize*500
				ax_nomerge.scatter(y,x,s=1,c='g')

			rho=projection_functions.read_cube(files_nomerge[-i-1])
			flat=np.sum(rho,2)/len(rho[:,:,0])
			ax_nomerge.imshow(np.log10(flat),cmap='magma',vmin=vmin,vmax=vmax)

	plt.subplots_adjust(hspace=0,wspace=-0.9)
	axs=ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10,ax11,ax12
	for ax in axs:
		ax.set_yticks([])
		ax.set_xticks([])
		
	
