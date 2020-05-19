import matplotlib.pyplot as plt
import numpy as np
from pac import Plot

contact0N = np.loadtxt('/home/eg/codes/pacOutput/OTC-0N/5D50-contactTableRW.csv', delimiter=',')
contact1500N = np.loadtxt('/home/eg/codes/pacOutput/OTC-1500N/5D50-contactTableRW.csv', delimiter=',')

Plot.equalAreaProjection(contact0N)
Plot.equalAreaProjection(contact1500N)


