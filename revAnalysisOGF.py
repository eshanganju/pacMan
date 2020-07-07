# Importing files
from pac import Reader                      # Reads files, cuts into smaller pieces
from pac import Filter                      # Filters files using NLM filter
from pac import Segment                     # Binarizatio and WS
from pac import Measure                     # Calculates particle size, mprphology, contact, breakage
from pac import Plot                        # Plotting functions
import time                                 # The fourth dimension
import matplotlib.pyplot as plt             # matplotlib
import skimage.external.tifffile as tf      # scikit-image
import numpy as np                          # numpy

inputFolderLocation = '/home/eg/codes/pacInput/OGF-0N/'
outputFolderLocation = '/home/eg/codes/pacOutput/REV/OGF-0N/'
originalGSDLocation = '/home/eg/codes/pacInput/originalGSD/ogfOrig.csv' # Original GSD location

# Data details 0N:
dataName = 'ogf-0N'
measuredVoidRatioSample = 0.635                                     # Void ratio measured from 1D compression experiment
d50 = 0.62                                                          # D50 in mm - original gradation
cal = 0.01193                                                       # calibration from CT mm/voxel
zCenter = 508                                                       # Voxel units - center of slice
yCenter = 437                                                       # Voxel units - vertical center
xCenter = 490                                                       # Voxel units - horizontal center
originalGSD = np.loadtxt( originalGSDLocation , delimiter=',' )     # Original GSD

analyzeTotalVol = False
analyzeRevSizes = True

if analyzeTotalVol == True :
	edgeLength = 9*d50 # Max edge length in D50s

	gliMap = Reader.readTiffFileSequence( inputFolderLocation, zCenter, yCenter, xCenter, edgeLength, cal, invertImageData=False)
	binMap, edMap, edPeakMap, labMap = Segment.obtainLabelledMapUsingITKWS( gliMap , measuredVoidRatio=measuredVoidRatioSample , outputLocation=outputFolderLocation )
	correctedLabMap = Segment.fixErrorsInSegmentation( labMap , pad=2)
	noEdgeCorrectedLabMap = Segment.removeEdgeLabels( correctedLabMap )

	# Save files:
	gliName = outputFolderLocation + 'gliMap.tiff'
	binName = outputFolderLocation + 'binMap.tiff'
	labName = outputFolderLocation + 'labMap.tiff'
	correctedLabName = outputFolderLocation + 'corLabMap.tiff'
	noEdgeCorrectedLabName = outputFolderLocation + 'noEdgeCorLabMap.tiff'

	tf.imsave(gliName,gliMap.astype('uint32'))
	tf.imsave(binName,binMap.astype('uint32'))
	tf.imsave(labName,labMap.astype('uint32'))
	tf.imsave(correctedLabName , correctedLabMap.astype('uint32'))
	tf.imsave(noEdgeCorrectedLabName , noEdgeCorrectedLabMap.astype('uint32'))

if analyzeRevSizes == True :
	if analyzeTotalVol == False:
		
		binName = outputFolderLocation + 'binMap.tiff'
		binMap = tf.imread(binName).astype('uint32')

		gliName = outputFolderLocation + 'gliMap.tiff'
		gliMap = tf.imread(gliName).astype('uint32')
		
		
	sizeRange = np.arange(8, 9, 1)
	zCenterSu = binMap.shape[0]//2
	yCenterSu = binMap.shape[1]//2
	xCenterSu = binMap.shape[2]//2
	
	sizeList = []
	voidRatioList = []
	gsdList = []

	for i in sizeRange:
		length = i*d50
		sizeList.append(i)

		print('\nAnalyzing subregion of size : ' + str(i) + 'd50')

		upSlice = int(zCenterSu + round( ( length / 2 ) / cal ))
		lowSlice = int(zCenterSu - round( ( length / 2 ) / cal ))
		upRow = int(yCenterSu + round( ( length / 2 ) / cal ))
		lowRow = int(yCenterSu - round( ( length / 2 ) / cal ))
		upCol = int(xCenterSu + round( ( length / 2 ) / cal ))
		lowCol = int(xCenterSu - round( ( length / 2 ) / cal ))
		
		subRegionBinMap = binMap[lowSlice:upSlice,lowRow:upRow,lowCol:upCol]
		subVoidRatio = Segment.calcVoidRatio(subRegionBinMap)
		voidRatioList.append(subVoidRatio)

		subRegionGliMap = gliMap[lowSlice:upSlice,lowRow:upRow,lowCol:upCol]
		binMask, edMap, edPeaksMap, subRegionLabMap = Segment.obtainLabelledMapUsingITKWS( subRegionGliMap, knownThreshold = 7533,outputLocation=outputFolderLocation )
		subRegionCorrectedLabMap = Segment.fixErrorsInSegmentation( subRegionLabMap , pad=2)
		#subRegionNoEdgeCorrectedLabMap = Segment.removeEdgeLabels( subRegionCorrectedLabMap )
		#xx, xx, xx, gsd = Measure.gsd( subRegionNoEdgeCorrectedLabMap , calib=cal )
		#np.savetxt((outputFolderLocation+ str(np.round(i)) +'D50-gsdPCAmin.csv'), gsd, delimiter=',')

		contactTableRW, contactTableITK = Measure.contactNormalsSpam(subRegionCorrectedLabMap)
		np.savetxt((outputFolderLocation+ str(np.round(i)) +'D50-contactTableRW.csv'), contactTableRW, delimiter=',')    # Contact table RW
		np.savetxt((outputFolderLocation+ str(np.round(i)) +'D50-contactTableITK.csv'), contactTableITK, delimiter=',')  # Contact table ITK

	#plt.plot(sizeList,voidRatioList)
	#plt.ylim([0,2])
	#plt.xlim([0,10])
	#plt.show()
