#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 14:32:24 2019

@author: lewisprole
"""


import numpy as np
import matplotlib.pyplot as plt 
import struct
import binascii





'''arepo initial conditions writer

N-dimentional arrays must be in v=xv,yv,zv form'''






def header(sofar,npart,massarr,time,redshift,flag_sfr,flag_feedback,
           npartTotal,flag_cooling,num_files,boxsize,cos1,cos2,
           hubble_param,flag_stellarage,flag_metals,npartHighword,
           flag_entropy,flag_dp,flag_1pt,scalefactor):
    '''function to write the header'''
    
    len_bytes=256
    sofar=struct.pack('i',8)       #create the tag for the header 
    sofar=sofar+bytes('HEAD',encoding='utf-8')
    sofar=sofar+struct.pack('i',256+8)
    sofar=sofar+struct.pack('i',8)
    
    npart=struct.pack('iiiiii',*npart)    #convert header info into bytes 
    massarr=struct.pack('dddddd',*massarr)
    time=struct.pack('d',time)
    redshift=struct.pack('d',redshift)
    flag_sfr=struct.pack('i',flag_sfr)
    flag_feedback=struct.pack('i',flag_feedback)
    npartTotal=struct.pack('iiiiii',*npartTotal)
    flag_cooling=struct.pack('i',flag_cooling)
    num_files=struct.pack('i',num_files)
    boxsize=struct.pack('d',boxsize)
    cos1=struct.pack('d',cos1)
    cos2=struct.pack('d',cos2)
    hubble_param=struct.pack('d',hubble_param)
    flag_stellarage=struct.pack('i',flag_stellarage)
    flag_metals=struct.pack('i',flag_metals)
    npartHighword=struct.pack('iiiiii',*npartHighword)
    flag_entropy=struct.pack('i', flag_entropy)
    flag_dp=struct.pack('i',flag_dp)
    flag_1pt=struct.pack('i',flag_1pt)
    scalefactor=struct.pack('i',scalefactor)
    
    bytes_left =256 - 6*4 - 6*8 -2*8 - 10*4 -4*8 - 12*4 #account for extra junk
    num_left=int(bytes_left/4)
    dummy=np.zeros(num_left).astype(int)
    filler=struct.pack('i'*num_left,*dummy)
    
    sofar=sofar+struct.pack('i',256)  #create header block
    sofar=sofar+npart+massarr+time+redshift+flag_sfr+flag_feedback+npartTotal
    sofar=sofar+flag_cooling +num_files + boxsize + cos1 + cos2 +hubble_param
    sofar=sofar+flag_stellarage+flag_metals+npartHighword+flag_entropy
    sofar=sofar+flag_dp+flag_1pt+scalefactor+filler
    sofar=sofar+struct.pack('i',256)
    
    return sofar


#h=header([],(10,0,0,0,0,0),(5,0,0,0,0,0),6,7,8,9,(9,0,0,0,0,0),1,2,3,4,5,6,7,8,(1,2,3,4,5,6),1,2,3,4)




def tagger(sofar,data,name,dtype,no_axis):
    '''function to write the tag before a data block'''
    a=struct.pack('i',8)
    

    
    if dtype =='i':
        f=4
    if dtype=='d':
        f=8
        
    if no_axis>1:
        N=int(len(data[0])*no_axis*f)
    else:
        N=int(len(data)*f) #number of args in data * no.bytes per number 
    
    n=struct.pack('i',N+8) #convert to bytes (+8 for the 2 boarders)
    
    print('writing %s'%name)
    if len(sofar)>0:
        sofar=sofar+a
        
    else:
        sofar=a
        
    sofar=sofar+bytes(name,encoding='utf-8')

    sofar=sofar+n

    sofar=sofar+a

    
    return sofar




def blocker(sofar,data,dtype,num_axis):
    '''function to write the block of data'''
    
    if num_axis>1:
        N=len(data[0])*num_axis
        
        if dtype == 'i': #choose dtype for conversion to binary 
            fmt='i'*N
            no_bytes=N*4
            
        if dtype == 'd':
            fmt ='d'*N
            no_bytes=N*8
        
        data=np.stack(data,axis=1) #want 1D array in form (x1,y1,z1,x2,y2,z2...) 
        data=np.reshape(data,(-1,N))
        
        block=struct.pack(fmt,*data[0]) #convert to binary
    
    else:   
            
        N=len(data)
        
        if dtype == 'i': #choose dtype for conversion to binary 
            fmt='i'*N
            no_bytes=N*4
            
        if dtype == 'd':
            fmt ='d'*N
            no_bytes=N*8
        
        block=struct.pack(fmt,*data) #convert to binary
    
    
    n=struct.pack('i',no_bytes) 
    print(no_bytes)
    
    sofar=sofar+n   #writing block 
    sofar=sofar+block
    sofar=sofar+n
    
    return sofar





def ghost(sofar,dtype,npart,name,num_axis):
    '''creates empty data arrays when data is absent'''
    zeros=np.zeros((num_axis,npart))
    sofar=tagger(sofar,zeros,name,dtype)
    sofar=blocker(sofar,zeros,dtype)
    return sofar 
        





def tag_block(sofar,data,name,dtype,no_axis):
    '''writes tag and data block of one data set into binary'''
    sofar=tagger(sofar,data,name,dtype,no_axis)
    sofar=blocker(sofar,data,dtype,no_axis)
    return sofar





#test
sofar,npart,massarr,time,redshift,flag_sfr,flag_feedback,npartTotal,flag_cooling,num_files,boxsize,cos1,cos2,hubble_param,flag_stellarage,flag_metals,npartHighword,flag_entropy,flag_dp,flag_1pt,scalefactor=[],(10,0,0,0,0,0),(5,0,0,0,0,0),6,7,8,9,(9,0,0,0,0,0),1,2,3,4,5,6,7,8,(1,2,3,4,5,6),1,2,3,4

sofar=header(sofar,npart,massarr,time,redshift,flag_sfr,flag_feedback,
           npartTotal,flag_cooling,num_files,boxsize,cos1,cos2,
           hubble_param,flag_stellarage,flag_metals,npartHighword,
           flag_entropy,flag_dp,flag_1pt,scalefactor)
sofar=tag_block(sofar,v,'VEL ','d',3)
sofar=tag_block(sofar,x,'POS ','d',3)
sofar=tag_block(sofar,ids,'ID  ','i',1)


    