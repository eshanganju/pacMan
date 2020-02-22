import spam.label as slab
import skimage.external.tifffile as tf

labMap = tf.imread('/home/eg/codes/pacOutput/specialSPAM/otc-0n-labelledMapITKWatershed-NE.tiff').astype(int)

cLabMap = slab.detectAndFixOversegmentation(labMap,nVoxThreshold=25)

tf.imsave('/home/eg/codes/pacOutput/specialSPAM/otc-0n-labelledMapITKWatershed-NE-Corrected.tiff',cLabMap)