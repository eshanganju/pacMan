# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np

numPts=1000

a = (np.random.rand(numPts, 3).reshape(numPts,3))
a_mean = np.ones_like(a)

meanX = np.average(a[:,0])
meanY = np.average(a[:,1])
meanZ = np.average(a[:,2])

a_mean[:,0] = a_mean[:,0] * meanX
a_mean[:,1] = a_mean[:,1] * meanY
a_mean[:,2] = a_mean[:,2] * meanZ
a_c = a-a_mean

I_a = np.zeros((3,3))
for i in range(0,numPts):
    I_a[0,0] = I_a[0,0] + ( a_c[i,1]**2 + a_c[i,2]**2 )
    I_a[1,1] = I_a[1,1] + ( a_c[i,0]**2 + a_c[i,2]**2 )
    I_a[2,2] = I_a[2,2] + ( a_c[i,0]**2 + a_c[i,1]**2 )
    
    I_a[0,1] = I_a[0,1] + ( a_c[i,0] * a_c[i,1] )
    I_a[0,2] = I_a[0,2] + ( a_c[i,0] * a_c[i,2] )
    I_a[1,2] = I_a[1,2] + ( a_c[i,1] * a_c[i,2] )
    
    I_a[1,0] = I_a[0,1]
    I_a[2,0] = I_a[0,2]
    I_a[2,1] = I_a[1,2]
    
I_a[0,1] = -I_a[0,1]
I_a[0,2] = -I_a[0,2]
I_a[1,2] = -I_a[1,2]

I_a[1,0] = I_a[0,1]
I_a[2,0] = I_a[0,2]
I_a[2,1] = I_a[1,2]

COVMat_a = np.cov( a_c.T )

EVal_I_a, EVec_I_a = np.linalg.eig( I_a )
EVal_COV_a, EVec_COV_a = np.linalg.eig( COVMat_a )

print('a')
print(a)
print('\n')

print('a_c')
print(a_c)
print('\n')

print('I_a')
print(I_a)
print('\n')

print('COVMat_a')
print(COVMat_a)
print('\n')

print('EVec_I_a')
print(EVec_I_a)
print('\n')

print('EVec_COV_a')
print(EVec_COV_a)
print('\n')


# Rotate particle using I_a
rot_a_I_a =  np.matmul( EVec_I_a.T, a_c.T )
lenX_I_a = rot_a_I_a[0,:].max() - rot_a_I_a[0,:].min()
lenY_I_a = rot_a_I_a[1,:].max() - rot_a_I_a[1,:].min()
lenZ_I_a = rot_a_I_a[2,:].max() - rot_a_I_a[2,:].min()
len_I_a = np.array( [ lenX_I_a, lenY_I_a, lenZ_I_a ] ).reshape(1,3)
print('Lengths I_a')
print(len_I_a)

rot_a_COV_a =  np.matmul( EVec_COV_a.T, a_c.T )
lenX_COV_a = rot_a_COV_a[0,:].max() - rot_a_COV_a[0,:].min()
lenY_COV_a = rot_a_COV_a[1,:].max() - rot_a_COV_a[1,:].min()
lenZ_COV_a = rot_a_COV_a[2,:].max() - rot_a_COV_a[2,:].min()
len_COV_a = np.array( [ lenX_COV_a, lenY_COV_a, lenZ_COV_a ] ).reshape(1,3)
print('Lengths COV_a')
print(len_COV_a)

