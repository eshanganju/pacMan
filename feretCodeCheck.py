import tifffile as tiffy
from pac import Measure
import numpy as np

scanList = ['OTC_25_mid',
            'OTC_25_tip',
            'OTC_25_top',
            '2QR_25_mid',
            '2QR_25_tip',
            '2QR_25_top']
folderLoc = '/home/eg/codes/pacOutput/cone/'

for scan in scanList:
    for num in range(0,10):
        fileLoc = folderLoc + scan + '/' + scan + '-' + str(int(num)) + '-noEdgeCorrectedLabelMap.tif'
        labMap = tiffy.imread(fileLoc).astype('uint16')
        particleSizeSummary = np.zeros((labMap.max(),8))
        particleSizeSummary = Measure.getParticleSize(labMap, calibrationFactor=1,
                                                      sampleName=scan+'-'+str(int(num))+'-feretIncl',
                                                      saveData=True, outputDir=folderLoc+scan+'/')

