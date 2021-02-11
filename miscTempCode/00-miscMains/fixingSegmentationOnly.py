import skimage.external.tifffile as tf
import classes.Segment as seg



labMap = tf.imread('/home/eg/codes/pacOutput/specialSPAM/otc-0n-labelledMapITKWatershed-small.tiff').astype(int)

# labMap = labMap[0:201, 0:201, 0:201]

# tf.imsave('dataRandom.tiff', labMap)

print('Number of labels before correction: %d' % labMap.max())

clm = seg.fixErrorsInSegmentation(labMap)

print('\nNumber of labels after correction: %d' % clm.max())

tf.imsave('/home/eg/codes/pacOutput/specialSPAM/otc-0n-small-labMapITKWatershed-corrected.tiff',clm)
