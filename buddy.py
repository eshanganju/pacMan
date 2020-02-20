from classes import Segment
from classes import Measure
import skimage.external.tifffile as tiffy
import matplotlib.pyplot as plt

outputLocation = '/home/eg/codes/pacOutput/'
fileName ='otc0N'

currentVoidRatio = 0.5405

s = Segment.Segment()

gli = tiffy.imread('/home/eg/codes/pacOutput/otc0N.tiff')
gliMax = (2**16)-1

fltr = tiffy.imread('/home/eg/codes/pacOutput/otc0N-filtered.tiff')
#fltr = Filter.filterDenoiseNlm(fileName,gli,gliMax,outputLocation)

binMap, threshold = s.binarizeAccordingToDensity( fltr, currentVoidRatio, outputLocation, fileName )
mapLocation = outputLocation + fileName + '-binMap.tiff'
tiffy.imsave( mapLocation, binMap)

edMap = s.obtainEuclidDistanceMap(binMap)
mapLocation = outputLocation + fileName + '-edMap.tiff'
tiffy.imsave( mapLocation, edMap)

edPeaksMarkerMap = s.obtainLocalMaximaMarkers(edMap)
mapLocation = outputLocation + fileName + '-edPeaksMap.tiff'
tiffy.imsave( mapLocation, edPeaksMarkerMap)

labelledMap = s.obtainLabelledMapUsingWaterShedAlgorithm( binMap,  edMap,  edPeaksMarkerMap )
mapLocation = outputLocation + fileName + '-labelledMap.tiff'
tiffy.imsave( mapLocation, labelledMap )

#cLabelledMap = s.fixErrorsInSegmentation(llabelledMap)
#gsdTable = m.measureParticleSizeDistribution(cLabelledMap, outputLocation)
#contactTable = m.measureContactNormalsSpam(binMap, clabelledMap, outputLocation)

#plt.figure()
#plt.imshow( binMap[binMap.shape[0] // 2 ], cmap ='gray' )

#plt.figure()
#plt.imshow( edMap[edMap.shape[0] // 2 ], cmap ='gray' )

#plt.figure()
#plt.imshow( edPeaksMarkerMap[edPeaksMarkerMap.shape[0] // 2 ], cmap ='gray' )

#plt.figure()
#plt.imshow( labelledMap[labelledMap.shape[0] // 2 ], cmap ='rainbow' )

#plt.show()
