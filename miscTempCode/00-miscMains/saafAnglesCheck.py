# -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 15:35:21 2021

@author: LENOVO USER
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

N = 1000

phi = np.zeros((N,1))
theta = np.zeros((N,1))

for i in range(0,N):
    n=i+1
    h = -1 + 2*(n-1)/(N-1) 
    theta[i] = np.arccos(h)

for i in range(1,N-1):
    n=i+1
    h = -1 + 2*(n-1)/(N-1) 
    val = phi[ i-1 ] + 3.6/( N )**0.5 * ( 1 / ( 1-h**2 )**0.5 )
    phi[i] = val%(2*np.pi)
        
z = np.cos(phi)
y = np.sin(phi) * np.sin(theta)
x = np.sin(phi) * np.cos(theta)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x, y, z) 