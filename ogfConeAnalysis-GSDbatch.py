import numpy as np
import matplotlib.pyplot as plt
import glob
import skimage.external.tifffile as tf

from classes import Measure


inputOutputDataLocation = '/home/eg/codes/pacInput/zaManualSegmentation2020-03-04/editedLabelFiles/'

readFromFile = False

if readFromFile == False:
    searchString = inputOutputDataLocation + '/*tif'
    fileList = glob.glob(searchString)
    print(len(fileList))
    print(fileList)

else:
    fileList = ['/home/eg/codes/pacInput/zaManualSegmentation2020-03-04/editedLabelFiles/ogf-25kPa-tip-spamLabelledMap-edited.tif']

limitFileNumber=1000000
numFilesToRead = min(limitFileNumber,len(fileList))
print('Yahoo')
for fileNum in range(0,numFilesToRead):
    labelledData = tf.imread(fileList[fileNum])
    particleSizeData = Measure.measureParticleSize(labelledData)
    fileName = (fileList[fileNum].split('/')[-1]).split('.')[0]
    completeFileName = inputOutputDataLocation + fileName + '.csv'
    np.savetxt(completeFileName,particleSizeData,delimiter=',')


print('--dobby--')

