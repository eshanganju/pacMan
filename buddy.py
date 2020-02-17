from classes import Segment
import skimage.external.tifffile as tiffy

fltr = tiffy.imread('C:/Users/eganj/gitHub/pacOutput/otc-filtered.tiff')
outputLocation = 'C:/Users/eganj/gitHub/pacOutput/'
fileName = 'otc'

s = Segment.Segment()

# s.binarizeAccordingToOtsu(fltr,outputFile,name) 
s.binarizeAccordingToUserThreshold(fltr,outputLocation,fileName) 

