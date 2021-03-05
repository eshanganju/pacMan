# -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 17:00:46 2021

@author: LENOVO USER
"""
import numpy as np
from numpy import arange, pi, sin, cos, arccos
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

n = 400

goldenRatio = (1 + 5**0.5)/2

i = arange(0, n)

theta = 2* pi * i / goldenRatio

phi = arccos(1 - 2*(i+0.5)/n)

x, y, z = cos(theta) * sin(phi), sin(theta) * sin(phi), cos(phi)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x, y, z) 

a = np.zeros((x.shape[0],3))

a[:,0] = x
a[:,1] = y
a[:,2] = z

name='fiboSpiral-' + str(n) + '.csv'

np.savetxt(name,a,delimiter=',')