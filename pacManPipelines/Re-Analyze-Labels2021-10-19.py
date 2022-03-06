""" 2021-09-03 @EG
Main code for the individual segmentation of particles - FGM datasets

This code uses the Segment module to convert the binary data to labelled map

The input needs to be in tiff format, preferably 8 bit
The output will be in tiff as well - labelled map

The particles will be analyzed

"""

#Imports
import tifffile as tf
from pac import Measure
from pac import Segment
import numpy as np

# Output
ofl = '/home/chawlahpc2adm/Desktop/Daniel/RAMpowder/Output/'
dataName = 'Al_643x643x969_zyx_5h_rr09'

# Read file
labelMap = tf.imread('/home/chawlahpc2adm/Desktop/Daniel/RAMpowder/Input/Al_643x643x969_zyx_5h_rr09-noEdgeCorrectedLabelMap_4.tif')

# Find missing labels
incorrectLabelsArray = np.arange(1,labelMap.max()+1,dtype='uint16') 
correctLabelsArray = np.unique(labelMap)

missingLabels = np.zeros((incorrectLabelsArray.shape[0] - correctLabelsArray.shape[0] + 1 ),dtype='uint16')
missingLabelPos = 0

for num in incorrectLabelsArray:
	if num in correctLabelsArray:None
	
	else:
		print(str(num) + ' missing')
		missingLabels[missingLabelPos]=num
		missingLabelPos = missingLabelPos+1

missingLabels = np.flip( missingLabels )

for missingLabel in missingLabels: labelMap = Segment.removeLabelAndUpdate(labelMap,missingLabel)
tf.imsave('Al_643x643x969_zyx_5h_rr09_UPDATE.tif',labelMap)
	
Measure.getParticleSizeArray( labelledMapForParticleSizeAnalysis = labelMap, 
								calibrationFactor=1,	# mm per voxel 
								saveData=True, 
								sampleName=dataName, 
								outputDir=ofl)
								

		


	                                     

