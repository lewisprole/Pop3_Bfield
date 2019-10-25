#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 12:46:38 2019

@author: lewisprole
"""

import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits import mplot3d
from matplotlib import colors
import code_units
import astropy.constants as ap 


def zero_vel(n):
    '''function to give all cells 0 velocity'''
    vx=np.zeros(n)
    vy=np.zeros(n)
    vz=np.zeros(n)
    v=vx,vy,vz
    return v


def rot_sphere(size,x,M,r,B):
    '''function to give sphere solid body rotation,
    all parameters to be given in cgs'''
    
    x,y,z=x[0],x[1],x[2]
    
    #if KE_rotation = A * U_grav then w=root(BM/r^3)
    
    G=ap.G.cgs.value
    w=np.sqrt(B*2*G*M)/r**3 
    mid=size/2
    
    distx=(mid-x) #distances from z axis of rotation 
    disty=(mid-y)
    dist=np.sqrt((mid-x)**2+(mid-y)**2)


    v_rot=w*dist #rotational velocity perpendicular to r 
    
    theta=np.arctan(disty/distx) #rosolve x and y velocities 
    

    mask_ypos=np.where(disty>0)
    
    mask_xpos=np.where(distx>0)
    mask_both_pos=np.intersect1d(mask_ypos,mask_xpos)
    
    mask_yneg=np.where(disty<0)
    mask_xneg=np.where(distx<0)
    mask_both_neg=np.intersect1d(mask_yneg,mask_xneg)
    
    ypos_xneg=np.intersect1d(mask_ypos,mask_xneg)
    xpos_yneg=np.intersect1d(mask_xpos,mask_yneg)
    
    vx=v_rot*np.cos(theta)
    vy=v_rot*np.sin(theta)
    vz=0
    
    
    vx[mask_both_pos]=-vx[mask_both_pos]
    vy[mask_both_neg]=-vy[mask_both_neg]
    vy[ypos_xneg]=-vy[ypos_xneg]
    vx[xpos_yneg]=-vx[xpos_yneg]
    
    
    rs = np.sqrt((mid-x)**2+(mid-y)**2+(mid-z)**2)  #no rotation outside cloud 
    mask=np.where(rs>r)
    vx[mask]=0 
    vy[mask]=0
    vy[mask]=0
    
    

    
    return vx,vy,vz
    


def vary_rotation(size,x,B,m):
    '''different angular momentum depending on mass within your radius,
    all parameters to be given in cgs'''
    x,y,z=x[0],x[1],x[2]
    mid=size/2
    distx=(mid-x) #distances from z axis of rotation 
    disty=(mid-y)
    dist=np.sqrt((mid-x)**2+(mid-y)**2)
    rs=np.sqrt((mid-x)**2+(mid-y)**2+(mid-z)**2)
    
    inM=0
    E=np.zeros_like(x)
    v_rotation=np.zeros_like(x)
    args=dist.argsort()
    
    
    
    G=ap.G.cgs.value
    for i in args:
        inM+=m[i]
        if rs[i]>0:
            E[i]=G*inM*m[i]/rs[i]
        else:
            E[i]=0
        if dist[i]>0:
            w=np.sqrt(B*2*E[i]/(inM*dist[i]**2))
        else:
            w=0
        v_rotation[i]=w*dist[i]

    
    theta=np.arctan(disty/distx) #rosolve x and y velocities 
    theta[np.where(distx==0)]=np.pi/2
        
    
    mask_ypos=np.where(disty>0)  #different rules for different combinations of x/y
    mask_xpos=np.where(distx>0)
    mask_both_pos=np.intersect1d(mask_ypos,mask_xpos) #both positive 
    
    mask_yneg=np.where(disty<0)
    mask_xneg=np.where(distx<0)
    mask_both_neg=np.intersect1d(mask_yneg,mask_xneg) #both negative 
    
    ypos_xneg=np.intersect1d(mask_ypos,mask_xneg) #mix of pos/neg
    xpos_yneg=np.intersect1d(mask_xpos,mask_yneg)        
    
    x0=np.where(distx==0)
    x0ypos=np.intersect1d(mask_ypos,x0)
    x0yneg=np.intersect1d(mask_yneg,x0)
    
    
    vx=v_rotation*np.sin(theta) #resolve velocit
    vy=v_rotation*np.cos(theta)
    vz=np.zeros_like(vx)
    
    
    vx[mask_both_pos]=-vx[mask_both_pos]
    vy[mask_both_neg]=-vy[mask_both_neg]
    vy[ypos_xneg]=-vy[ypos_xneg]
    vx[xpos_yneg]=-vx[xpos_yneg]
    vx[x0ypos]=v_rotation[x0ypos]
    vx[x0yneg]=-v_rotation[x0yneg]
    vy[x0]=0
    
    return vx,vy,vz
        
        


#v=zero_vel(10000)


#def max_boltz(n,T,v):
#    '''function to give Maxwell-Boltzmann velocity distribution'''
#    kb=1.38e-23
#    
#    dn=4*np.pi*n*(m/(2*np.pi*kb*T))**(1/2)*v**2*np.exp((-m*v**2)/2*kb*T)
#    
#    return dn
#
#def mb_dist():
    
    
