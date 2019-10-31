#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 12:21:26 2019

@author: lewisprole
"""

import numpy as np
import astropy.constants as ap 

def BE_radius(method,M_or_rho,T_ism):
    '''calculate radius of bonner ebert sphere given either
the mass of the sphere or the ambient background density, also
returns either rho_out or M_cloud (whichever was not given)'''
    G=ap.G.cgs.value
    kb=ap.k_B.cgs.value
    mp=ap.m_p.to('g').value
    G=ap.G.cgs.value
    mu=2.4
    c_s=np.sqrt(kb*T_ism/(mu*mp))
    if method=='mass':
        R=M_or_rho*G/(2.4*c_s**2)
        ungiven=(1.18*c_s**3/(G**(3/2)*M_or_rho))**2 #rho
    if method=='rho':
        R=0.49*c_s/(np.sqrt(G*M_or_rho))
        ungiven=1.18*c_s**3/(G**(3/2)*M_or_rho**(1/2)) #M
    return R,ungiven 

def temp(M,rho):
    '''calculate temperature of ISM from sphere mass and outter density'''
    G=ap.G.cgs.value
    kb=ap.k_B.cgs.value
    mp=ap.m_p.to('g').value
    G=ap.G.cgs.value
    mu=2.4
    c_s=(M*G**(3/2)*rho**(1/2)/1.18)**(1/3)
    T=c_s**2*mu*mp/kb
    return T


def checkR(R,T_ism):
    '''function to check calculations'''
    G=ap.G.cgs.value
    kb=ap.k_B.cgs.value
    mp=ap.m_p.to('g').value
    G=ap.G.cgs.value
    mu=2.4
    c_s=np.sqrt(kb*T_ism/(mu*mp))
    
    M=2.4*c_s**2*R/G
    rho=(0.49*c_s/(G**(1/2)*R))**2
    return M,rho    
