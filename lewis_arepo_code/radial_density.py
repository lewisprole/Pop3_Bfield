#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 11:07:49 2019

@author: lewisprole
"""

import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits import mplot3d
from matplotlib import colors

'''radial density distribution'''

def rho_dist(r,A,a,f,c ):
    '''density distribution within cloud as function of position r'''
    rho=A*(a*r)**f +c 
    return rho


def rhos(xs,ys,zs,x_size,y_size,z_size,R_cloud,A,a,f,c):
    '''set up desnity inside and outside cloud'''
    
    midx,midy,midz=int(x_size/2),int(y_size/2),int(z_size/2)
    
    rs=np.sqrt((midx-xs)**2+(midy-ys)**2+(midz-zs)**2)
    print(rs)
    
    rho=rho_dist(rs,A,a,f,c) #give all cells radial distribution
    
    mask=np.where(rs>R_cloud) #set outside sphere as constant
    bg=max(rho[mask])
    rho[mask]=bg
    
    return rho, rs

#rho,rs=rhos(x[0],x[1],x[2],1000,1000,1000,100,100,100**2,-2,0)
    

#fig = plt.figure()
#ax = plt.axes(projection='3d')
#im = ax.scatter3D(x,y,z, s=1,c=rho,cmap='hot')
#fig.colorbar(im)
#
#plt.plot(rs,rho,'x')

