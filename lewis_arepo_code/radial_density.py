#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 11:07:49 2019

@author: lewisprole
"""

import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits import mplot3d
from matplotlib import colors
import astropy.constants as ap
from scipy import interpolate
import math 
'''radial density distribution'''

def rho_dist(r,A,a,f,c ):
    '''density distribution within cloud as function of position r'''
    rho=A*(a*r)**f +c 
    return rho


def rhos(xs,ys,zs,x_size,y_size,z_size,R_cloud,A,a,f,c):
    '''set up desnity inside and outside cloud'''
    
    midx,midy,midz=int(x_size/2),int(y_size/2),int(z_size/2)
    
    rs=np.sqrt((midx-xs)**2+(midy-ys)**2+(midz-zs)**2)
    print(rs)
    
    rho=rho_dist(rs,A,a,f,c) #give all cells radial distribution
    
    mask=np.where(rs>R_cloud) #set outside sphere as constant
    bg=max(rho[mask])
    rho[mask]=bg
    
    return rho, rs

def BE_profile(xs,ys,zs,size,T,rho_bg):
    '''Bonnor Ebert sphere density profile'''
    mid=size/2
    kb=ap.k_B.cgs.value
    mp=ap.m_p.to('g').value
    G=ap.G.cgs.value
    mu=2.4

    c_s=np.sqrt(kb*T/(mu*mp)) #in cgs
    
    rs=np.sqrt((mid-xs)**2+(mid-ys)**2+(mid-zs)**2)
    mask=np.where(rs>0)
    rho=np.zeros_like(xs)
    rho[mask]=c_s**2/(2*np.pi*G*rs[mask]**2)
    
    mask=np.where(rs==0)
    rho[mask]=rho.max()
    
    R=0.49 *c_s/(G**(1/2)*rho_bg**(1/2))
    mask=np.where(rs>R)
    rho[mask]=rho_bg

    return rho

def non_crit_BE(x,y,z,size,T,n0,n_bg,R):
	mid=size/2
	kb=ap.k_B.cgs.value
	mp=ap.m_p.to('g').value
	G=ap.G.cgs.value
	mu=2.4
	c_s=np.sqrt(kb*T/(mu*mp))
	RS=np.sqrt((mid-x)**2+(mid-y)**2+(mid-z)**2)
	



	r_crit=1.2e4*ap.au.cgs.value
	rho_crit=3*c_s**2/(2*np.pi*G*r_crit**2)
	
	a=np.sqrt(c_s**2/(4*np.pi*G*n0))
	rho=n0/(1+(RS**2/(3*a**2)))
	mask=np.where(RS>R)
	rho[mask]=n_bg
	return rho,RS
	

	#phi=0
	#rs=np.linspace(0.1,R,1e10)	
	#args=rs.argsort()
	#Min=4/3*np.pi*rs[args[0]]**3*n0
	#rhos=np.zeros_like(RS)
	#rhos[0]=n0
	#for i in range (len(rs)-1):
	#	I=args[i]
	#	J=args[i+1]
	#	rho=n0*np.exp(-G*Min/(rs[I]*c_s**2))
	#	rhos[i+1]=rho
	#	Min+=(rs[J]-rs[I]) *rho *4*np.pi*rs[I]**2 
	#mask=np.where(RS>R)
	#rhos[mask]=n_bg
	#return rhos


	#rs=np.linspace(0.1,R,1e5)
	#z=rs/c_s *np.sqrt(4*np.pi*G*n0)
	#dz=z[1]-z[0]
	#w=1/6 * z**2 - 1/(5*math.factorial(4))*z**4 + 8/(21*math.factorial(6))*z**6 
	#dw=np.array([])
	#for i in range (len(w)-1):
	#	dw=np.append(dw,w[i+1]-w[i])
	#dw=np.append(dw,0)
	#y=dw/dz	
	#y=z**2*y
	#dy=np.array([])
	#for i in range (len(y)-1):
	#	dy=np.append(dy,y[i+1]-y[i])
	#dy=np.append(dy,0)
	#rho=np.zeros_like(RS)
	#mask=np.where(RS<=R)
	#rho[mask]=n0/z[mask]**2 * dy[mask]/dz
	#mask=np.where(RS>R)
	#rho[mask]=n_bg 	
	#return rho



	#set up BE spline 
	#rs=np.linspace(1e-10,R,1e8)
	#inM=n0*4/3 *np.pi*rs[0]**3
	#rho=np.zeros_like(rs)
	#rho[0]=n0
	#rho_i=n0
	#for i in range(len(rs)-1):
	#	dr=(rs[i+1]-rs[i])

	#	k1=-G*inM*rho_i/(c_s**2*rs[i]**2) *dr
	#	
	#	

	#	r=rs[i]+dr/2
	#	rho2=rho_i+k1/2
	#	dm2=inM +4*np.pi*rs[i]**2 *rho2 *(dr/2)
	#	k2=-G*dm2*rho2/(c_s**2*r**2) *dr

	#	rho3=rho_i+k2/2
	#	dm3=inM +4*np.pi*rs[i]**2 *rho3 *(dr/2)
	#	k3=-G*dm2*rho3/(c_s**2*r**2) *dr

	#	r=rs[i+1]
	#	rho4=rho_i+k3
	#	dm4=inM +4*np.pi*rs[i]**2 *rho4 *(dr)
	#	k4=-G*dm4*rho4/(c_s**2*r**2) *dr

	#	rho_i+=1/6 * (k1 +2*k2 +2*k3 * k4)

	#	rho[i+1]=rho_i
	#	inM += 4*np.pi*rs[i]**2 *rho_i *dr
	#return rho 

	#	rho_i+=drho
	#	rho[i+1]=rho_i
	#	inM+=4*np.pi*rs[i]**2*rho_i*(rs[i+1]-rs[i])
	#spine=interpolate.interp1d(rs, rho, kind='cubic')
	#rhos=np.ones_like(RS)*n_bg
	#mask=np.where(RS<R)
	#rhos[mask]=spine(RS[mask])
	#return rhos,spine
	

	#args=rs.argsort()
	#inM=n0*4/3*np.pi*rs[args[0]]**3 
	#rho_i=n0
	#rho=np.array([n0])
	#for i in range(len(args)-1):
	#	I=args[i]
	#	J=args[i+1]
	#	drho=-G*inM*rho_i*(rs[J]-rs[I])/(c_s**2*rs[I]**2)
	#	rho_i+=drho
	#	rho=np.append(rho,rho_i)
	#	inM+= 4*np.pi*rs[I]**2 * rho_i *(rs[J]-rs[I])
	#
	#mask=np.where(rs>R)
	#rho[mask]=n_bg
	#return rho

def tmass(r,A,a,f,c):
    '''total mass of the cloud'''
    R=np.linspace(0,r,100)
    dR=R[1]-R[0]
    n=0
    M=0
    R_current=R[1]
    while R_current<r:
        dm=rho_dist(R_current,A,a,f,c)*2*np.pi*dR
        M+=dm
        n+=1
        R_current=R[n]
        
    return M
        
        
#rho,rs=rhos(x[0],x[1],x[2],1000,1000,1000,100,100,100**2,-2,0)
    

#fig = plt.figure()
#ax = plt.axes(projection='3d')
#im = ax.scatter3D(x,y,z, s=1,c=rho,cmap='hot')
#fig.colorbar(im)
#
#plt.plot(rs,rho,'x')

