#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 14:20:02 2019

@author: lewisprole
"""

import numpy as np
import code_units
import astropy.constants as ap 

'''mass'''

def bonnor_ebert(size,x,v,v_out,T,r):
    '''bonner ebert radial density profile, balance gravity/pressure'''
    
    mid=int(size/2)
    x,y,z=x[0],x[1],x[2]
    rs=np.sqrt((mid-x)**2+(mid-y)**2+(mid-z)**2) #distances from center 
    
    kb=ap.k_B.cgs.value
    mp=ap.m_p.to('g').value
    G=ap.G.cgs.value
    mu=2.4

    c_s=np.sqrt(kb*T/(mu*mp)) #in cgs 
    print(c_s)
    
    rho=c_s**2/(2*np.pi*G*rs**2) #rs should already be in cm
    
    rho[np.isinf(rho)]=rho[~np.isinf(rho)].max()
    
    rho[np.where(rs>r)]=rho[np.where(rs<r)].min() #flat profile outside sphere 
    
    

    
    m=rho*v #in g  
    m[np.where(rs<r)]=rho[np.where(rs<r)]*v_out
    
    
    return m,rs,rho



 
