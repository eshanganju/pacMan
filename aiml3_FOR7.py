"""Extraction of LOF
"""

import numpy as np
import tifffile as tf

# Output:
ofl='/home/eg/Desktop/Al-AM-samples/output/finalAnalysis-7/'

# Limits
solidLimit = 8
voidLimit = 10

# Read EDMs, and corrected label maps
solidLabelMap = tf.imread('/home/eg/Desktop/Al-AM-samples/output/finalAnalysis-7/sample7-CROP-BIN-Solid161dist-watershed.tif').astype('uint16')
voidLabelMap = tf.imread('/home/eg/Desktop/Al-AM-samples/output/finalAnalysis-7/sample7-CROP-BIN-VOID160dist-watershed.tif').astype('uint16')

solidEDM = tf.imread('/home/eg/Desktop/Al-AM-samples/output/finalAnalysis-7/sample7-CROP-BIN-Solid161-dist.tif')
voidEDM = tf.imread('/home/eg/Desktop/Al-AM-samples/output/finalAnalysis-7/sample7-CROP-BIN-VOID160-dist.tif')


# Create LOF maps
lofVoid = np.zeros_like(voidLabelMap)
lofSolid = np.zeros_like(solidLabelMap)

# Analyzed label size for each void label
for voidLabel in range(1, voidLabelMap.max()+1):
	print("VOID")
	print('\tChecking void ' + str(voidLabel) + '/' +  str(voidLabelMap.max()))
	maxInscSph = voidEDM[np.where(voidLabelMap == voidLabel)].max()
	
	if maxInscSph <= voidLimit:
		lofVoid[np.where(voidLabelMap == voidLabel)] = 1

tf.imwrite( (ofl+'lof-void_' + str(voidLimit) + '.tif'), lofVoid)


# Analyzed label size for each solid label
for solidLabel in range(1, solidLabelMap.max()+1):
	print("SOLID")
	print('\tChecking solid ' + str(solidLabel) + '/' +  str(solidLabelMap.max()))
	maxInscSph = solidEDM[np.where(solidLabelMap == solidLabel)].max()
	
	if maxInscSph <= solidLimit:
		lofSolid[np.where(solidLabelMap == solidLabel)] = 1

tf.imwrite( (ofl+'lof-solid_' + str(solidLimit) + '.tif'), lofSolid)

combinedMap = lofVoid +  lofSolid
if combinedMap.max() > 1:
	print('Overlap')
	combinedMap[np.where(combinedMap>1)]=1

tf.imwrite((ofl+'lof-solid_' + str(solidLimit) +'-void_' + str(voidLimit) + '.tif'), combinedMap)