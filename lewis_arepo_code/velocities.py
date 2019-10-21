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

def zero_vel(n):
    '''function to give all cells 0 velocity'''
    vx=np.zeros(n)
    vy=np.zeros(n)
    vz=np.zeros(n)
    v=vx,vy,vz
    return v


def rot_sphere(size,x,M,r,B):
    '''function to give sphere solid body rotation'''
    x,y,z=x[0],x[1],x[2]
    
    #if KE_rotation = A * U_grav then w=root(BM/r^3)
    
    w=np.sqrt(B*M)/r**3 
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
    '''different angular momentum depending on mass within your radius'''
    x,y,z=x[0],x[1],x[2]
    mid=size/2
    distx=(mid-x) #distances from z axis of rotation 
    disty=(mid-y)
    dist=np.sqrt((mid-x)**2+(mid-y)**2)
    
    inM=0
    E=np.zeros_like(x)
    v_rotation=np.zeros_like(x)
    args=dist.argsort()
    
    
    
    
    for i in args:
        inM+=m[i]
        E[i]=inM*m[i]/dist[i]
        w= np.sqrt(B*2*E[i]/(inM*dist[i]**2))
        v_rotation[i]=w*dist[i]
    
    theta=np.arctan(disty/distx) #rosolve x and y velocities 
    
    mask_ypos=np.where(disty>0)  #different rules for different combinations of x/y
    mask_xpos=np.where(distx>0)
    mask_both_pos=np.intersect1d(mask_ypos,mask_xpos) #both positive 
    
    mask_yneg=np.where(disty<0)
    mask_xneg=np.where(distx<0)
    mask_both_neg=np.intersect1d(mask_yneg,mask_xneg) #both negative 
    
    ypos_xneg=np.intersect1d(mask_ypos,mask_xneg) #mix of pos/neg
    xpos_yneg=np.intersect1d(mask_xpos,mask_yneg)        
    vx=v_rotation*np.cos(theta) #resolve velocit
    vy=v_rotation*np.sin(theta)
    vz=0
    
    
    vx[mask_both_pos]=-vx[mask_both_pos]
    vy[mask_both_neg]=-vy[mask_both_neg]
    vy[ypos_xneg]=-vy[ypos_xneg]
    vx[xpos_yneg]=-vx[xpos_yneg]
    
    return vx,vy,vz,v_rotation
        
        


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
    
    
