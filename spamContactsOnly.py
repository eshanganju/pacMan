import numpy as np
import spam.plotting as splt
import spam.label as slab
import skimage.external.tifffile as tf
from matplotlib import cm
from classes import Reader

labMap = tf.imread('/home/eg/codes/pacOutput/specialSPAM/otc-1500n-labMapITKWatershed-NE.tiff').astype(int)
fgliMap = tf.imread('/home/eg/codes/pacOutput/otc-1500n-fgliMap.tiff').astype(int)

contVol, z, contTable, contLab = slab.labelledContacts(labMap) 
contOrt = slab.contactOrientationsAssembly(labMap,fgliMap,contLab,watershed='RW', NumberOfThreads=4)

# contOrt = np.genfromtxt('/home/eg/codes/pacOutput/specialSPAM/otc-500n-contOrt-RW.csv',delimiter=',')

contOrtOnly = contOrt[:,2:5]
tempContOrtOnly = np.zeros_like(contOrtOnly)
count = 0

for row in range(0,contOrtOnly.shape[0]):
	length = (contOrtOnly[row,0]**2+contOrtOnly[row,1]**2+contOrtOnly[row,2]**2)**0.5
	
	if length == 1:
		tempContOrtOnly[count] = contOrtOnly[row]
		count = count + 1

	elif length > 1:
		contOrtOnly[row,:] = contOrtOnly[row,:]/length
		tempContOrtOnly[count] = contOrtOnly[row]
		count = count + 1
		
		print('Length > 1')
		print(contOrtOnly[row])
		print(length)
		print('\n')
	
	elif length < 1:
		
		if length ==0:
			print('Length = 0')
			print(contOrtOnly[row])
			print(length)
			print('\n')
		
		else: 
			contOrtOnly[row,:] = contOrtOnly[row,:]/length
			print('Length < 1')
			print(contOrtOnly[row])
			print(length)
			print('\n')
			tempContOrtOnly[count] = contOrtOnly[row]
			count = count + 1

cleanContOrtOnly = tempContOrtOnly[0:count,:]

splt.plotOrientations(cleanContOrtOnly, binValueMax = 8, numberOfRings = 18, cmap = cm.Reds)

f1,f2,d = slab.fabricTensor(cleanContOrtOnly)
np.savetxt('/home/eg/codes/pacOutput/specialSPAM/otc-1500n-contOrt-RW.csv',contOrt,delimiter=',')
np.savetxt('/home/eg/codes/pacOutput/specialSPAM/otc-1500n-f1-RW.csv',f1,delimiter=',')
np.savetxt('/home/eg/codes/pacOutput/specialSPAM/otc-1500n-f2-RW.csv',f2,delimiter=',')


