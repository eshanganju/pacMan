import numpy as np
from uncertainties import unumpy as unp
from uncertainties import ufloat
from uncertainties.umath import *

fileName = '/home/eg/codes/pacOutput/OGF-1500N-8/6.0D50-contactTableRW.csv'

data = np.loadtxt(fileName, delimiter=',')

vectors = data[:,2:5]
uncertVectors = 0.26*(np.ones_like(vectors))

uncertVectorArray = unp.uarray(vectors,uncertVectors)

N11 = ufloat(0,0)
N12 = ufloat(0,0)
N13 = ufloat(0,0)
N21 = ufloat(0,0)
N22 = ufloat(0,0)
N23 = ufloat(0,0)
N31 = ufloat(0,0)
N32 = ufloat(0,0)
N33 = ufloat(0,0)

for i in range(0,uncertVectorArray.shape[0]):
    N11 = N11 + (uncertVectorArray[i,0])*(uncertVectorArray[i,0])
    N12 = N12 + (uncertVectorArray[i,0])*(uncertVectorArray[i,1])
    N13 = N13 + (uncertVectorArray[i,0])*(uncertVectorArray[i,2])
    N21 = N21 + (uncertVectorArray[i,1])*(uncertVectorArray[i,0])
    N22 = N22 + (uncertVectorArray[i,1])*(uncertVectorArray[i,1])
    N23 = N23 + (uncertVectorArray[i,1])*(uncertVectorArray[i,2])
    N31 = N31 + (uncertVectorArray[i,2])*(uncertVectorArray[i,0])
    N32 = N32 + (uncertVectorArray[i,2])*(uncertVectorArray[i,1])
    N33 = N33 + (uncertVectorArray[i,2])*(uncertVectorArray[i,2])

N11 = N11 / uncertVectorArray.shape[0]
N12 = N12 / uncertVectorArray.shape[0]
N13 = N13 / uncertVectorArray.shape[0]
N21 = N21 / uncertVectorArray.shape[0]
N22 = N22 / uncertVectorArray.shape[0]
N23 = N23 / uncertVectorArray.shape[0]
N31 = N31 / uncertVectorArray.shape[0]
N32 = N32 / uncertVectorArray.shape[0]
N33 = N33 / uncertVectorArray.shape[0]

traceF = N11+N22+N33

F11 = 15/2 * ( N11 - (1/3)*traceF )
F12 = 15/2 * ( N12 )
F13 = 15/2 * ( N13 )
F21 = 15/2 * ( N21 )
F22 = 15/2 * ( N22 - (1/3)*traceF )
F23 = 15/2 * ( N23 )
F31 = 15/2 * ( N31 )
F32 = 15/2 * ( N32 )
F33 = 15/2 * ( N33 - (1/3)*traceF )

Fq = ((3/2)*( F11*F11 + F12*F12 + F13*F13 + F21*F21 + F22*F22 + F23*F23 + F31*F31 + F32*F32 + F33*F33))**0.5

print('Fq = ' + str(Fq))
