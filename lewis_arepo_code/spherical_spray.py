#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 17:10:43 2019

@author: lewisprole
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

'''script to set up initial positions'''

def sphere_fill(n,r,x_size,y_size,z_size):
    '''spray in n particles within sphere radius r'''
    midx,midy,midz=int(x_size/2),int(y_size/2),int(z_size/2)
    xs,ys=np.random.randint(midx-r,midx+r,n),np.random.randint(midy-r,midy+r,n)
    zs=np.random.randint(midz-r,midz+r,n)
    rs=np.sqrt((midx-xs)**2+(midy-ys)**2+(midz-zs)**2)
    
    mask=np.where(rs<r)
    xs=xs[mask]
    ys=ys[mask]
    zs=zs[mask]

    
    return xs, ys, zs
    


def bg_fill(n_bg,r,x_size,y_size,z_size):
    '''sprays in background particles'''
    midx,midy,midz=int(x_size/2),int(y_size/2),int(z_size/2)
    xs,ys, =np.random.randint(0,x_size,n_bg),np.random.randint(0,y_size,n_bg)
    zs=np.random.randint(0,z_size,n_bg)
    rs=np.sqrt((midx-xs)**2+(midy-ys)**2+(midz-zs)**2)
    mask=np.where(rs>r)
    xs=xs[mask]
    ys=ys[mask]
    zs=zs[mask]
    return xs, ys, zs


def spherical_cloud(n,n_bg,r,x_size,y_size,z_size):
    xs_c,ys_c,zs_c=sphere_fill(n,r,x_size,y_size,z_size)
    xs_bg,ys_bg,zs_bg=bg_fill(n_bg,r,x_size,y_size,z_size)
    
    return np.append(xs_c,xs_bg),np.append(ys_c,ys_bg),np.append(zs_c,zs_bg)




#
def plotter(n,n_bg,r,x_size,y_size,z_size):
#    xs_c,ys_c,zs_c=sphere_fill(n,r,x_size,y_size,z_size)
#    xs_bg,ys_bg,zs_bg=bg_fill(n,r,x_size,y_size,z_size)
    x,y,z=spherical_cloud(n,n_bg,r,x_size,y_size,z_size)
    
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    
#    ax.scatter3D(xs_c,ys_c,zs_c, s=1,cmap='Greens')
#    ax.scatter3D(xs_bg,ys_bg,zs_bg, s=1,cmap='Reds')
    ax.scatter3D(x,y,z, s=1,cmap='Reds')



#plotter(1000,1000,100,1000,1000,1000)
x,y,z=spherical_cloud(1000,1000,100,1000,1000,1000)
x=x,y,z
ids=np.linspace(1,1001,1001).astype(int)
