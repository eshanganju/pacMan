import matplotlib.pyplot as plt
import numpy as np
from pac import Plot
from pac import Measure

for i in range(1,4):
    if i == 1:
        print('\n\nFor OTC:')
        print('---------------*')
        contact0N = np.loadtxt('/home/eg/codes/pacOutput/OTC-0N/5D50-contactTableRW.csv', delimiter=',')
        contact500N = np.loadtxt('/home/eg/codes/pacOutput/OTC-500N/5D50-contactTableRW.csv', delimiter=',')
        contact1500N = np.loadtxt('/home/eg/codes/pacOutput/OTC-1500N/5D50-contactTableRW.csv', delimiter=',')

    elif i == 2:
        print('\n\nFor OGF:')
        print('---------------*')
        contact0N = np.loadtxt('/home/eg/codes/pacOutput/OGF-0N/5D50-contactTableRW.csv', delimiter=',')
        contact500N = np.loadtxt('/home/eg/codes/pacOutput/OGF-500N/5D50-contactTableRW.csv', delimiter=',')
        contact1500N = np.loadtxt('/home/eg/codes/pacOutput/OGF-1500N/5D50-contactTableRW.csv', delimiter=',')

    elif i == 3:
        print('\n\nFor 2QR:')
        print('---------------*')
        contact0N = np.loadtxt('/home/eg/codes/pacOutput/2QR-0N/5D50-contactTableRW.csv', delimiter=',')
        contact500N = np.loadtxt('/home/eg/codes/pacOutput/2QR-500N/5D50-contactTableRW.csv', delimiter=',')
        contact1500N = np.loadtxt('/home/eg/codes/pacOutput/2QR-1500N/5D50-contactTableRW.csv', delimiter=',')

    #Plot.equalAreaProjection(contact0N,nr=9)
    #Plot.equalAreaProjection(contact500N,nr=9)
    #Plot.equalAreaProjection(contact1500N,nr=9)

    N0N , F0N , Fq0N = Measure.fabricVariables( contact0N )
    N500N , F500N , Fq500N = Measure.fabricVariables( contact500N )
    N1500N , F1500N , Fq1500N = Measure.fabricVariables( contact1500N )


    print('N0N')
    print(N0N)
    print('F0N')
    print(F0N)
    print('Fq0N')
    print(Fq0N)
    print('\n')

    print('N500N')
    print(N500N)
    print('F500N')
    print(F500N)
    print('Fq500N')
    print(Fq500N)
    print('\n')

    print('N1500N')
    print(N1500N)
    print('F1500N')
    print(F1500N)
    print('Fq1500N')
    print(Fq1500N)
