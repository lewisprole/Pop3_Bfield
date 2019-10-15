#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 10:24:37 2019

@author: lewisprole
"""

import numpy as np
import matplotlib.pyplot as plt
import struct
import binascii
from  matplotlib import colors

#filename='/Users/lewisprole/Documents/PhD/Pop3_Bfield/snapshot_004'





filename='snapshot_3636'
filename='snapshot_344'
#filename="arepo_input.txt"

with open(filename, mode='rb') as file:
    data = file.read()



def load_head(Data,tag_start):
    '''function to read the tag before the data block'''
    
    a=Data[tag_start:tag_start+4] #initial tag boarder is 4 bytes (1 int32)
    
    name=Data[tag_start+4:tag_start+8] #4 bytes for the 4 letter name 
    
    name=str(name,'utf-8') #convert bytes into string
    
    length=Data[tag_start+8:tag_start+12] #length of data block given by an int32 (4bytes)
    if len(length)>0:
        stopper=1
        length_unpack=struct.unpack('i',length)[0] #convert from binary 
        length_real=length_unpack-8 #data length: tag length minus boarders
        a_check=Data[tag_start+12:tag_start+16] #final tag boarder is 4 bytes (1 int32)
        
    else:
        stopper=0
        length_real=0
        a_check=0
        name, length_real=0,0
    
    
    return name, length_real, tag_start+16,stopper









def load_data(Data,length_real,block_start,name):
    '''function for reading the data block after the tag'''

    a=Data[block_start:block_start+4] #boarderis 4 bytes = 1 int32
    
    block=Data[block_start+4:block_start+length_real+4]  #length of data in bytes given by load_head
    
    a_check=Data[block_start+length_real+4:block_start+length_real+4+4]
    print(struct.unpack('i',a),struct.unpack('i',a_check))
    if a==a_check:
        print(name+': matching i/f')
    else:
        print(name +': non-matching i/f')

    return block, block_start+length_real+4+4








def reader(Data):
    '''function to loop through tag/block pairs and pass info into class 'a' '''
    class a(): #empty class to give all the info to
        pass
    
    
    tag_start=0 #start at the beginning of the snapshot file 

    stopper=1 #used for ending the loop
    while stopper>0:
        
        name, length_real, block_start, stopper=load_head(Data,tag_start)
        
        if stopper>0:
            block, tag_start=load_data(Data,length_real,block_start,name)
            
            
    
            if name=='HEAD':
        
                s=0
                npart=struct.unpack('iiiiii',block[s:s+24])
                s+=24
                massarr=struct.unpack('dddddd',block[s:s+48])
                s+=48
                time=struct.unpack('d',block[s:s+8])
                s+=8
                redshift=struct.unpack('d',block[s:s+8])
                s+=8
                flag_sfr=struct.unpack('i',block[s:s+4])
                s+=4
                flag_feedback=struct.unpack('i',block[s:s+4])
                s+=4
                npartTotal=struct.unpack('iiiiii',block[s:s+24])
                s+=24
                flag_cooling=struct.unpack('i',block[s:s+4])
                s+=4
                num_files=struct.unpack('i',block[s:s+4])
                s+=4
                boxsize=struct.unpack('d',block[s:s+8])
                s+=8
                cos1=struct.unpack('d',block[s:s+8])
                s+=8
                cos2=struct.unpack('d',block[s:s+8])
                s+=8
                hubble_param=struct.unpack('d',block[s:s+8])
                s+=8
                flag_stellarage=struct.unpack('i',block[s:s+4])
                s+=4
                flag_metals=struct.unpack('i',block[s:s+4])
                s+=4
                npartHighword=struct.unpack('iiiiii',block[s:s+24])
                s+=24
                flag_entropy=struct.unpack('i',block[s:s+4])
                s+=4
                flag_dp=struct.unpack('i',block[s:s+4])
                s+=4
                flag_1pt=struct.unpack('i',block[s:s+4])
                s+=4
                scalefactor=struct.unpack('i',block[s:s+4])
                s+=4
                
                
                num_left= 256 - 24-48-8-8-4-4-24-4-4-8-8-8-8-4-4-24-4-4-4-4    
                
                filler=struct.unpack('i'*int(num_left/4), block[s:s+num_left])
                
                
                
                print("npart array:", npart)
                print("mass array:", massarr)
                print("time in codeunits:", time)
                print("Total npart:", npartTotal)
                print("Header finished")

                N = int(sum(npart)) - npart[4] # type 4 is reserved for TRACER_MC
                ngas = npart[0]
                nsink = npart[5]

                #
                # Initialise our struct that holds the data                
                
                a.npart = npart
                a.ngas = ngas
                a.nsink = nsink
                a.N = N
                a.time = time
                a.redshift=redshift
                a.flag_sfr=a.flag_sfrflag_sfr
                a.flag_feedback=flag_feedback
                a.npartTotal=npartTotal
                a.flag_cooling=flag_cooling
                a.num_files=num_files
                a.boxsize=boxsize
                a.cos1=cos1
                a.cos2=cos2
                a.hubble_param=hubble_param
                a.flag_stellarage=flag_stellarage
                a.flag_metals=flag_metals
                a.npartHighword=npartHighword
                a.flag_entropy=flag_entropy
                a.flag_dp=flag_dp
                a.flag_1pt=flag_1pt
                a.scalefactor=scalefactor
                
                
                
                #
                # Set the units here
                a.unit_leng_cm = 1.0e+17
                a.unit_mass_g = 1.991e33
                a.unit_time_s = 2.7436898e12 
                
                

            
        if name=='POS ':
            print("Reading positions")
            n=int(length_real/8) #vel data stored as ints (4bytes per number)
            fmt='d'*n
            pos=struct.unpack(fmt,block)
            pos=np.reshape(pos,(-1,3))
            a.x=pos[:,0]
            a.y=pos[:,1]
            a.z=pos[:,2]
        
        elif name=='VEL ':
            print("Reading velocities")
            n=int(length_real/8)
            fmt='d'*n
            vel=struct.unpack(fmt,block)
            vel=np.reshape(vel,(-1,3))
            a.vx=vel[:,0]
            a.vy=vel[:,1]
            a.vz=vel[:,2]
            
        elif name=='ID  ':
            print("Reading IDs")
            n=int(length_real/4)
            fmt='i'*n
            ids=struct.unpack(fmt,block)
            a.ids=ids
            
        elif name=='MASS':
            print("Reading masses")
            n=int(length_real/8)
            fmt='d'*n
            mass=struct.unpack(fmt,block)
            a.mass=mass
            
        elif name=='U   ':
            print("Reading U")
            n=int(length_real/8)
            fmt='d'*n
            u=struct.unpack(fmt,block)
            a.u=u        
            
        elif name=='RHO ':
            print("Reading densities")
            n=int(length_real/8)
            fmt='d'*n
            rho=struct.unpack(fmt,block)
            a.rho=rho        
            
        elif name=='POT ':
            print("Reading potentials")
            n=int(length_real/8)
            fmt='d'*n
            pot=struct.unpack(fmt,block)
            a.pot=pot    
            
        elif name=='ACCE':
            print("Reading accelerations")
            n=int(length_real/8)
            fmt='d'*n
            acce=struct.unpack(fmt,block)
            a.acce = np.reshape(acce,(-1,3))

        
        elif name=='TSTP':
            print("Reading TSTP")
            n=int(length_real/8)
            fmt='d'*n
            tstp=struct.unpack(fmt,block)
            a.tstp=tstp
            
        elif name=='DIVV':
            print("Reading DIVV")
            n=int(length_real/8)
            fmt='d'*n
            divv=struct.unpack(fmt,block)
            a.divv=divv
            
        elif name=='SOFT':
            print("Reading SOFT")
            n=int(length_real/8)
            fmt='d'*n
            soft=struct.unpack(fmt,block)
            a.soft=soft
            
        elif name=='PEAK':
            print("Reading the peaks of the potential")
            n=int(length_real/4)
            fmt='i'*n
            peak=struct.unpack(fmt,block)
            a.peak=peak

                
            
    return a


def hexer(a):
    plt.hexbin(a.y,a.x,gridsize=1000,Norm=colors.LogNorm())
#    plt.xlim(2.1,2.6)
#    plt.ylim(2.1,2.6)
    plt.colorbar()
    
    
#a=reader(data)
m=reader(data)



