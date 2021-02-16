"""Analysis of  large amount of crushing
"""

from pac import Reader
from pac import Segment
from pac import Measure
from pac import Plot

v = 0

if v==0:

if v==1:

if v==2:

gliMap = Reader.extractSubregionFromTiffSequence()
binMap = Segment.binarizeAccordingToOtsu()
edtMap = Segment.obtainEuclidDistanceMap()
edpMap = Segment.obtainLocalMaximaMarkers()
labMap = Segment.segmentUsingWatershed()
corLabMap = Segment.fixErrorsInSegmentation()
fixCorLabMap = Segment.removeSmallParticles()
noEdgeFixCorLabMap = Segment.removeEdgeLabels()
psArray = Measure.getParticleSizeArray()
psdFeret = Measure.getParticleSizeDistribution()
relBrkHardin = Measure.getRelativeBreakageHardin()

