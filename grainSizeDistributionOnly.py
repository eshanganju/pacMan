import skimage.external.tifffile as tf
from classes import Measure
import numpy as np

m = Measure.Measure()

labMap = tf.imread('/home/eg/codes/pacOutput/specialSPAM/otc-1500n-labMapITKWatershed-NE.tiff')
psds = m.measureParticleSizeDistribution(labMap)

np.savetxt('/home/eg/codes/pacOutput/specialSPAM/otc-1500n-gsd.csv',psds,delimiter=',')