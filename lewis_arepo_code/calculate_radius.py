#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 12:21:26 2019

@author: lewisprole
"""

import numpy as np
import astropy.constants as ap 

def BE_radius(M,T_ism):
    '''calculate radius of bonner ebert sphere'''
    G=ap.G.cgs.value
    kb=ap.k_B.cgs.value
    mp=ap.m_p.to('g').value
    G=ap.G.cgs.value
    mu=2.4
    c_s=np.sqrt(kb*T_ism/(mu*mp))
    
    R=M*G/(2.4*c_s**2)
    return R 
    