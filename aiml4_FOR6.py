"""Variation of LOF and porosity
"""

import tifffile as tf
import numpy as np

sampleName = 'sample6'

ifl = '/home/eg/Desktop/Al-AM-samples/output/finalAnalysis-6/'
ofl = ifl

lofBinFile = '/home/eg/Desktop/Al-AM-samples/output/finalAnalysis-6/lof-solid_6-void_10_BIN-ERODE1_DILATE1.tif'
porosityBinFile = '/home/eg/Desktop/Al-AM-samples/output/finalAnalysis-6/sample6-CROP-MEDIAN_1-Void_OTSU.tif'

lofBin = tf.imread(lofBinFile)
poreBin = tf.imread(porosityBinFile)

lofBin = ( lofBin // lofBin.max() )
poreBin = ( poreBin // poreBin.max() )

numSlices = lofBin.shape[0]
dataPercentLofAndPorosity = np.zeros((numSlices,3)).astype('float')

sliceVolume = float(1.0 * lofBin.shape[1] * lofBin.shape[2])

for sliceNumber in range(0,numSlices):
	print('Analyzing slice ' + str(sliceNumber) + '/' + str(numSlices-1))
	lofVolume = lofBin[sliceNumber].sum()
	poreVolume = poreBin[sliceNumber].sum()
	
	percentLof = lofVolume/sliceVolume 
	percentPore = poreVolume/sliceVolume

	dataPercentLofAndPorosity[ sliceNumber, 0 ] = sliceNumber
	dataPercentLofAndPorosity[ sliceNumber, 1 ] = percentLof
	dataPercentLofAndPorosity[ sliceNumber, 2 ] = percentPore

np.savetxt( ( ofl + sampleName + '-VolumeData.csv' ), dataPercentLofAndPorosity, delimiter=',')
