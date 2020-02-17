from classes import Segment
import skimage.external.tifffile as tiffy

fltr = tiffy.imread('C:/Users/eganj/gitHub/pacOutput/otc-filtered.tiff')
outputLocation = 'C:/Users/eganj/gitHub/pacOutput/'
fileName = 'otc'

currentVoidRatio = 0.5404

s = Segment.Segment()

# binMap = s.binarizeAccordingToOtsu(fltr,outputFile,name) 
# binMap = s.binarizeAccordingToUserThreshold(fltr, outputLocation, fileName) 
binMap = s.binarizeAccordingToDensity(fltr, currentVoidRatio, outputLocation, fileName)
