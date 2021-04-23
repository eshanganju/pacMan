'''
Plotting lambert projects using PAC
written by @eg 2020-08-27
'''

import numpy as np
from pac import Plot
import matplotlib.pyplot as plt

fileLocation = '/home/eg/codes/pacOutput\
                /2QR-500N-2/6.0D50-contactTableRW.csv'

contTable = np.loadtxt(fileLocation, delimiter=',')

Plot.equalAreaProjection(contTable,
                         plotToPresentOption='both',
                         xSize=9,
                         ySize=4)
