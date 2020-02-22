import numpy as np
from classes import Reader
from classes import Filter
import skimage.external.tifffile as tf

r = Reader.Reader()
f = Filter.Filter()

currentVoidRatio = 0.5170
cntrZ = (20+1006)//2
cntrY = 461
cntrX = 502
lngt = 5.5
calib = 0.01193
bitDepth = 16
gliMax = 2**bitDepth-1

folderLocation = '/home/eg/codes/pacInput/OTC-1500N'
outputLocation = '/home/eg/codes/pacOutput/'
fileName = 'otc-1500n'

gliMap = r.readTiffFileSequence(folderLocation,cntrZ,cntrY,cntrX,lngt,calib)
gliMapName = outputLocation + fileName + '-gliMap.tiff'
tf.imsave(gliMapName,gliMap)


fgliMap = f.filterDenoiseNlm(fileName,gliMap,gliMax,outputLocation)
fgliMapName = outputLocation + fileName + '-fgliMap.tiff'
tf.imsave(fgliMapName,fgliMap)
