"""Pipeline for analysis of GB distance distribution for EBSD-DCT comparison
"""

import numpy as np
import tifffile as tf
from pac import Segment

# Global
ofl = '/home/eg/Desktop/GrainBoundaryOverlapAnalysis/Comparison-AlignedAndCropped/'


# Orient GBs


# Read-in oriented GB maps and binarize to 0 and 1
referenceMap = tf.imread('/home/eg/Desktop/GrainBoundaryOverlapAnalysis/Comparison-AlignedAndCropped/SEM-GB-Invert.tif')//255
gbMap_HDCT = tf.imread('/home/eg/Desktop/GrainBoundaryOverlapAnalysis/Comparison-AlignedAndCropped/HDCT-GB.tif')//255
gbMap_CDCT2 = tf.imread('/home/eg/Desktop/GrainBoundaryOverlapAnalysis/Comparison-AlignedAndCropped/CDCT2-GB.tif')//255

# Compute ED map of the reference map
edmReferenceMap = Segment.obtainEuclidDistanceMap( binaryMapForEDM=referenceMap, scaleUp = int(1), saveImg=True, 
													sampleName='SEM-Reference', 
													outputDir= ofl)

# Get multiplication of gb maps with ref maps
multi_refMap_HDCT = (edmReferenceMap * gbMap_HDCT).astype('float')
multi_refMap_CDCT2 = (edmReferenceMap * gbMap_CDCT2).astype('float')

# Save multiplication of gb maps with ref maps
tf.imwrite(ofl+'SEM-EDT_HDCT_product.tif',multi_refMap_HDCT.astype('uint8'))
tf.imwrite(ofl+'SEM-EDT_CDCT2_product.tif',multi_refMap_CDCT2.astype('uint8'))

# Get distance values for plots
multi_refMap_HDCT[np.where(multi_refMap_HDCT==0)]=np.nan
multi_refMap_CDCT2[np.where(multi_refMap_CDCT2==0)]=np.nan

flat_multi_refMap_HDCT = np.ndarray.flatten(multi_refMap_HDCT)
flat_multi_refMap_CDCT2 = np.ndarray.flatten(multi_refMap_CDCT2)

flat_multi_refMap_HDCT = flat_multi_refMap_HDCT[~np.isnan(flat_multi_refMap_HDCT)]
flat_multi_refMap_CDCT2 = flat_multi_refMap_CDCT2[~np.isnan(flat_multi_refMap_CDCT2)]

np.savetxt(ofl+'HDCT_GB_DIST.csv',flat_multi_refMap_HDCT,delimiter=',')
np.savetxt(ofl+'CDCT2_GB_DIST.csv',flat_multi_refMap_CDCT2,delimiter=',')