'''
This code is to test functionality of Filter using NLM
'''

from pac import Filter
import tifffile as tf

gliMap = tf.imread('scan0_5.tif')

bitDepthCurrent = 16
patchSizeCurrent = 5
patchDistanceCurrent = 15

outputLoc = ''
currentSampleName = 'scan0_5'
					

filteredGLIMap = Filter.filterUsingNlm( gli=gliMap,
										bitDepth=16,
										pSize=5,
										pDistance=7,
										saveImg=True,
										outputDir=outputLoc,
										sampleName=currentSampleName )