"""Pipeline To Compute selective coordination numbers
"""

from pac import Measure
import tifffile as tf

# These need to be of the same size
dataset1 =
dataset2 =

incrementMatrix = np.zeros_like(dataset2)
incrementMatrix[np.where(dataset2!=0)]=dataset1.max()
dataset2 = dataset2 + incrementMatrix

combinedDataset = dataset1 + dataset2

currentSampleName='trial'
ofl=''

Measure._getSelectiveCoordinationNumberList( labelMap=combinedDataset,
											 labelBoundary=dataset1.max()+1,
											 outputDir=ofl,sampleName=currentSampleName,
											 saveData=True )
