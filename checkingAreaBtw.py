import numpy as np
from pac import Measure

gsdC = np.genfromtxt('/home/eg/codes/pacInput/originalGSD/otcCurr.csv',delimiter=',')
gsdO = np.genfromtxt('/home/eg/codes/pacInput/originalGSD/otcOrig.csv',delimiter=',')

orig, curr, ult, Br = Measure.relativeBreakage(gsdC,gsdO,0.901,2.6)
