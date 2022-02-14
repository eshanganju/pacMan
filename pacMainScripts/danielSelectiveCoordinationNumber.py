"""Pipeline To Compute selective coordination numbers
"""

from pac import Measure
import tifffile as tf

# These need to be of the same size (ZYX)
dataset1 = tf.imread()
dataset2 = tf.imread()

incrementMatrix = np.zeros_like( dataset2 )
incrementMatrix[ np.where( dataset2!=0 ) ] = dataset1.max()
dataset2 = dataset2 + incrementMatrix

combinedDataset = dataset1 + dataset2

# TEST:
#combinedDataset=tf.imread('/home/chawlahpc2adm/Desktop/EG-NiCrC/pacMan/volumeAll.tif')

# sampleName
currentSampleName='Sample1'

# Output folder location
ofl=''

Carr = Measure._getSelectiveCoordinationNumberList( labelMap=combinedDataset,
											 labelBoundary=dataset1.max()-1,
											 outputDir=ofl,
											 sampleName=currentSampleName,
											 saveData=True )
											 
Carr					 
print(Carr)
											 
