from classes import Segment
import skimage.external.tifffile as tiffy
import matplotlib.pyplot as plt

fltr = tiffy.imread('C:/Users/eganj/gitHub/pacOutput/otc-filtered.tiff')
outputLocation = 'C:/Users/eganj/gitHub/pacOutput/'
fileName = 'otc'

currentVoidRatio = 0.5404

s = Segment.Segment()

binMap, threshold = s.binarizeAccordingToOtsu( fltr,
											   outputLocation,
											   fileName) 

binMap, threshold = s.binarizeAccordingToUserThreshold( fltr, 
														outputLocation, 
														fileName ) 

binMap, threshold = s.binarizeAccordingToDensity( fltr, 
												  currentVoidRatio, 
												  outputLocation, 
												  fileName)

edMap = s.obtainEuclidDistanceMap(binMap)
edPeaksMarkerMap = s.obtainLocalMaximaMarkers(edMap)

labelledMap = s.obtainLabelledMapUsingWaterShedAlgorithm( binMap, 
														  edMap, 
														  edPeaksMarkerMap )

plt.figure()
plt.imshow( binMap[binMap.shape[0] // 2 ], cmap ='gray' )

plt.figure()
plt.imshow( edMap[edMap.shape[0] // 2 ], cmap ='gray' )

plt.figure()
plt.imshow( edPeaksMarkerMap[edPeaksMarkerMap.shape[0] // 2 ], cmap ='gray' )

plt.figure()
plt.imshow( labelledMap[labelledMap.shape[0] // 2 ], cmap ='rainbow' )