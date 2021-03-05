# -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 16:47:52 2021

@author: LENOVO USER
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math

numOrt = 400
M = int(numOrt)*2
s = 3.6 / math.sqrt(M)
delta_z = 2 / float(M)
z = 1 - delta_z/2

longitude = 0

points = np.zeros( (numOrt,3) )

for k in range( numOrt ):
    r = math.sqrt( 1 - z*z )
    points[k,2] = math.cos( longitude ) * r     #X
    points[k,1] = math.sin( longitude ) * r     #Y
    points[k,0] = z                             #Z
    z = z - delta_z
    longitude = longitude + s/r
    

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(points[:,0], points[:,1], points[:,2]) 

name='saffAndK-' + str(numOrt) + '.csv'
np.savetxt(name,points,delimiter=',')
    
