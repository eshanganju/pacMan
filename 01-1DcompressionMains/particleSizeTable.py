'''
Objective:
    Compute the aspect ratio of the paticles.
    AR-1 - maxCentroidalAxisLength / medCentroidalAxesLength
    AR-2 - maxCentroidalAxisLength / minCentroidalAxesLength

Approach:
    Read noEdgeLabCorMao using tifffile
    Carry out particle size analysis using Measure
    Compute AR-1 and AR-2
    Save ration of AR1 and 2 as csv.
'''

import numpy as np
import math
from pac import Measure
import skimage.external.tifffile as tf


fileLocList = [ '/home/eg/codes/pacOutput/1D/2QR-0N/',
                '/home/eg/codes/pacOutput/1D/2QR-50N/',
                '/home/eg/codes/pacOutput/1D/2QR-100N/',
                '/home/eg/codes/pacOutput/1D/2QR-500N-2/',
                '/home/eg/codes/pacOutput/1D/2QR-1500N-2/',
                '/home/eg/codes/pacOutput/1D/OGF-0N/',
                '/home/eg/codes/pacOutput/1D/OGF-100N/',
                '/home/eg/codes/pacOutput/1D/OGF-500N/',
                '/home/eg/codes/pacOutput/1D/OGF-1500N-2/',
                '/home/eg/codes/pacOutput/1D/OTC-0N/',
                '/home/eg/codes/pacOutput/1D/OTC-500N/',
                '/home/eg/codes/pacOutput/1D/OTC-1500N/']


for fileLoc in fileLocList:
    fileName = fileLoc + 'noEdgeCorLabMap.tiff'
    noEdgeCorLabMap = tf.imread(fileName).astype('uint16')
    gss = Measure.getParticleSize(noEdgeCorLabMap, calibrationFactor = 0.01193)
    outputName = fileLoc + 'gss.csv'
    np.savetxt(outputName,gss,delimiter=',')


