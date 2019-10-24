#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 10:59:37 2019

@author: lewisprole
"""
from mpl_toolkits.mplot3d import Axes3D 
import matplotlib.pyplot as plt
import numpy as np


'''3D plotting'''

def plot3d(xs,ys,zs,size,c):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    if len(c)==len(xs):
        ax.scatter(xs, ys, zs, c=c,s=size)
    else:
        ax.scatter(xs, ys, zs, s=size)