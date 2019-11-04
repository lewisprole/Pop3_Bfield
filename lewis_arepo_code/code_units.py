#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 16:43:56 2019

@author: lewisprole
"""
import astropy.constants as ap 
import numpy as np

'''constants'''
M_cu=1.991e33 #g
d_cu =1e17 #cm
v_cu= 36447.268200 #cm/s
t_cu=d_cu/v_cu #s
rho_cu=M_cu / d_cu**3 #g/cm^3
B_cu=np.sqrt(rho_cu*v_cu**2)


kb=ap.k_B.cgs
mp=ap.m_p.to('g')