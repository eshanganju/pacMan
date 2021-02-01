import tifffile as tiffy
from pac import Measure
import numpy as np

fileLoc = '/home/eg/codes/pacOutput/cone/OTC_25_tip/OTC_25_tip-9-noEdgeCorrectedLabelMap.tif'
labMap = tiffy.imread(fileLoc).astype('uint16')

particleSizeSummary = np.zeros((labMap.max(),8))

particleSizeSummary = Measure.getParticleSize(labMap,
                                              calibrationFactor=0.01239099602,
                                              sampleName='OTC_25_tip',
                                              saveData=True)

