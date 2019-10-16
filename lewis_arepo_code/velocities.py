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
    
    
