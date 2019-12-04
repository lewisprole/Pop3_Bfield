#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 17:10:43 2019

@author: lewisprole
"""

import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits import mplot3d
import code_units
'''script to set up initial positions'''

def sphere_fill(n,r,x_size,y_size,z_size,precise):
    '''spray in n particles within sphere radius r'''
    midx,midy,midz=int(x_size/2),int(y_size/2),int(z_size/2)
    X=[]
    Y=[]
    Z=[]
    N=0
    if precise=='yes':
        while int(N)<n:
            xs=np.random.uniform(low=(midx-r),high=(midx+r),size=1)
            ys=np.random.uniform(low=(midy-r),high=(midy+r),size=1)
            zs=np.random.uniform(low=(midz-r),high=(midz+r),size=1)
            rs=np.sqrt((midx-xs)**2+(midy-ys)**2+(midz-zs)**2)
            if rs <= r:
                X=np.append(X,xs)
                Y=np.append(Y,ys)
                Z=np.append(Z,zs)
                N=len(X)
    else:
        xs=np.rand.uniform(low=(midx-r),high=(midx+r),size=n)
        ys=np.rand.uniform(low=(midx-r),high=(midx+r),size=n)       
        zs=np.rand.uniform(low=(midx-r),high=(midx+r),size=n)
        rs=np.sqrt((midx-xs)**2+(midy-ys)**2+(midz-zs)**2)
        mask=np.where(rs<=r)
        X=xs[mask]
        Y=ys[mask]
        Z=zs[mask]
    print('N_sphere: '+str(len(X)))
#    xs,ys=np.random.randint(midx-r,midx+r,size=n),np.random.randint(midy-r,midy+r,size=n)
#    zs=np.random.randint(midz-r,midz+r,size=n)
#    rs=np.sqrt((midx-xs)**2+(midy-ys)**2+(midz-zs)**2)
#    
#    mask=np.where(rs<r)
#    xs=xs[mask]
#    ys=ys[mask]
#    zs=zs[mask]
#
#    l=len(xs)
    return X,Y,Z
    


def bg_fill(n_bg,r,x_size,y_size,z_sizei,precise):
    '''sprays in background particles'''
    midx,midy,midz=int(x_size/2),int(y_size/2),int(z_size/2)
    X=[]
    Y=[]
    Z=[]
    N=0
    if precise=='yes':
        while int(N)<n_bg:
            xs=np.random.uniform(low=(0),high=(x_size),size=1)
            ys=np.random.uniform(low=(0),high=(y_size),size=1)
            zs=np.random.uniform(low=(0),high=(z_size),size=1)
            rs=np.sqrt((midx-xs)**2+(midy-ys)**2+(midz-zs)**2)
            if rs > r:
                X=np.append(X,xs)
                Y=np.append(Y,ys)
                Z=np.append(Z,zs)
                N=len(X)
    else:
        xs=np.rand.uniform(low=(midx-r),high=(midx+r),size=n)     
        ys=np.rand.uniform(low=(midx-r),high=(midx+r),size=n)
        zs=np.rand.uniform(low=(midx-r),high=(midx+r),size=n)
        rs=np.sqrt((midx-xs)**2+(midy-ys)**2+(midz-zs)**2)
        mask=np.where(rs>r)
        X=xs[mask]
        Y=ys[mask]
        Z=zs[mask]
    print('N_bg: '+str(len(X)))
#    xs,ys, =np.random.randint(0,x_size,size=n_bg),np.random.randint(0,y_size,size=n_bg)
#    zs=np.random.randint(0,z_size,size=n_bg)
#    rs=np.sqrt((midx-xs)**2+(midy-ys)**2+(midz-zs)**2)
#    mask=np.where(rs>r)
#    xs=xs[mask]
#    ys=ys[mask]
#    zs=zs[mask]
#    l=len(xs)
    return X,Y,Z


def spherical_cloud(n,n_bg,r,x_size,y_size,z_size,precise):
    xs_c,ys_c,zs_c,=sphere_fill(n,r,x_size,y_size,z_size,presice)
    xs_bg,ys_bg,zs_bg,=bg_fill(n_bg,r,x_size,y_size,z_size,precise)
    print('N_tot: '+str(len(xs_c)+len(xs_bg)))
    return np.append(xs_c,xs_bg),np.append(ys_c,ys_bg),np.append(zs_c,zs_bg)




def uniform_sphere(n,n_bg,r,size):
    '''n uniformly spaced positions on a grid within the sphere
    note n out can be different from the n given'''
    
    
    #fill in sphere
    mid=int(size/2)
    
    vol_cell=(4/3 * np.pi * r**3)/ n
    N=(size**3/vol_cell)
    N=round(N**(1/3))
    if N**3==n:
        N+=1  
    
    x=np.linspace(0,size,N)
    y=np.linspace(0,size,N)
    z=np.linspace(0,size,N)
    dx=x[1]-x[0]
    x, y, z = np.meshgrid(x, y, z)
    x, y, z = x.ravel(),y.ravel(),z.ravel()
    
    
    rs=np.sqrt((mid-x)**2 + (mid-y)**2 + (mid-z)**2)
    
    R=r 
    mask=np.where(rs<=R)
    x=x[mask]
    y=y[mask]
    z=z[mask]
    rs=rs[mask]
    
    
#    mask=np.array([])
#    if len(x)>n:
#        mask=[]
#        dif=len(x)-n
#        for i in range(dif):
#            mask=np.append(mask,rs.argmax())
#            rs[rs.argmax()]=0
#    
#    mask=mask.astype(int)
#    x=np.delete(x,mask)
#    y=np.delete(y,mask)
#    z=np.delete(z,mask)

    
    
    #now fill in background
    vol_cell_bg=(size**3 - (4/3 * np.pi * r**3 ))/n_bg
    N=round(size**3/vol_cell_bg)
    
    N=round(N**(1/3))
    if N**3==n_bg:
        N+=1
   
   
    buf=0.01*code_units.d_cu 
    xbg=np.linspace(buf,size-buf,N)
    ybg=np.linspace(buf,size-buf,N)
    zbg=np.linspace(buf,size-buf,N)
    xbg, ybg, zbg = np.meshgrid(xbg, ybg, zbg)
    xbg, ybg, zbg = xbg.ravel(),ybg.ravel(),zbg.ravel()
    

    rs=np.sqrt((mid-xbg)**2 + (mid-ybg)**2 + (mid-zbg)**2)
    
    R=r 
    mask=np.where(rs>R)
    xbg=xbg[mask]
    ybg=ybg[mask]
    zbg=zbg[mask]   
    rs=rs[mask]
    
    
    
#    mask=np.array([])
#    if len(xbg)>n_bg:
#        
#        dif=len(xbg)-n
#        
#        for i in range(dif):
#            
#            mask=np.append(mask,rs.argmax())
#            rs[rs.argmax()]=0
#    
#    mask=mask.astype(int)
#    xbg=np.delete(xbg,mask)
#    ybg=np.delete(ybg,mask)
#    zbg=np.delete(zbg,mask)
    
    
    
    #join sphere and background
    x=np.append(x,xbg)
    y=np.append(y,ybg)
    z=np.append(z,zbg)
    print(x.max()/code_units.d_cu)
    rs=np.sqrt((mid-x)**2 + (mid-y)**2 + (mid-z)**2)
    return x,y,z,vol_cell,vol_cell_bg


    


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
#x,y,z=spherical_cloud(5000,5000,100,1000,1000,1000)
#x=x,y,z
#ids=np.linspace(1,1001,1001).astype(int)
