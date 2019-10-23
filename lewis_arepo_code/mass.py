#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 14:20:02 2019

@author: lewisprole
"""

import numpy as np

'''mass'''

def bonnor_ebert(size,x,v,v_out,T,r):
    '''bonner ebert radial density profile, balance gravity/pressure'''
    
    mid=int(size/2)
    x,y,z=x[0],x[1],x[2]
    rs=np.sqrt((mid-x)**2+(mid-y)**2+(mid-z)**2) #distances from center 
    print(np.sort(rs))
    
    mu=2.4
    k_b=1.38e-23
    m_H=1.67e-27
    c_s=k_b*T/(mu*m_H)
    
    rho=c_s**2/(2*np.pi*rs**2)
    rho[np.isinf(rho)]=rho[~np.isinf(rho)].max()
    
    rho[np.where(rs>r)]=rho[np.where(rs<r)].min() #flat profile outside sphere 
    
    m=rho*v
    m[np.where(rs<r)]=rho[np.where(rs<r)]*v_out
    
    
    return m,rs,rho



 
