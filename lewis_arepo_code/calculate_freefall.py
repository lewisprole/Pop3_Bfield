import numpy as np
import astropy.constants as ap

def ff(m,x,r,size):
    '''free-fall time in code units'''
    #G=ap.G.cgs.value
    x,y,z=x[0],x[1],x[2]
    mid=size/2
    rs=np.sqrt((mid-x)**2+(mid-y)**2+(mid-z)**2)
    mask=np.asarray(np.where(rs<r)).astype(int)
    M=np.asarray(m)
    M=M[mask]
    M=sum(M[0])
    t=np.pi/2 * r**(3/2) / np.sqrt(2*M) #no G in sqrt because code units 
    return t 


