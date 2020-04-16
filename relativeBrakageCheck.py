import numpy as np
import matplotlib.pyplot as plt
from pac import Measure

'''
Import grain size distributions
Get relative breakage

'''

gsdOriginal = np.loadtxt('/home/eg/codes/pacInput/trialOriginalSort.csv',delimiter=',')
gsdCurrent  = np.loadtxt('/home/eg/codes/pacInput/trialCrushedSort.csv',delimiter=',')

gsdOrigNew,gsdCurNew,gsdUlt,Br = Measure.relativeBreakage(gsdOriginal,gsdCurrent,maxSize=0.4,fracDim=2.6)

np.savetxt('gsdOrigNew.csv',gsdOrigNew,delimiter=',')
np.savetxt('gsdCurNew.csv',gsdCurNew,delimiter=',')
np.savetxt('gsdUlt.csv',gsdUlt,delimiter=',')
