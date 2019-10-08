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

filename='/Users/lewisprole/Documents/PhD/Pop3_Bfield/snapshot_004'
data=np.loadtxt(filename)
#s = struct.Struct('<' + ' I 2s f')
#s.unpack(data)