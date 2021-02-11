'''
List of coordination numerbs for labelled map
    The image file is read from the corlabmap
    Measure gets the coordination number.
    The coordination number distribution is saved as text.
'''


import skimage.external.tifffile as tf
from pac import Measure
import numpy as np

'''
sample = '2QR-0N'
fileNameLoc = '/home/eg/codes/pacOutput/2QR-0N/corLabMap.tiff'
print('Checking for ' + sample)
corLabMap = tf.imread( fileNameLoc ).astype( 'uint32' )
coordinationNumberList = Measure.getCoordinationNumberList( corLabMap )
np.savetxt( (sample + 'cnList.txt'), coordinationNumberList, delimiter=',')

sample = '2QR-50N'
fileNameLoc = '/home/eg/codes/pacOutput/2QR-50N/corLabMap.tiff'
print('Checking for ' + sample)
corLabMap = tf.imread( fileNameLoc ).astype( 'uint32' )
coordinationNumberList = Measure.getCoordinationNumberList( corLabMap )
np.savetxt( (sample + 'cnList.txt'), coordinationNumberList, delimiter=',')

sample = '2QR-100N'
fileNameLoc = '/home/eg/codes/pacOutput/2QR-100N/corLabMap.tiff'
print('Checking for ' + sample)
corLabMap = tf.imread( fileNameLoc ).astype( 'uint32' )
coordinationNumberList = Measure.getCoordinationNumberList( corLabMap )
np.savetxt( (sample + 'cnList.txt'), coordinationNumberList, delimiter=',')

sample = '2QR-500N-2'
fileNameLoc = '/home/eg/codes/pacOutput/2QR-500N-2/corLabMap.tiff'
print('Checking for ' + sample)
corLabMap = tf.imread( fileNameLoc ).astype( 'uint32' )
coordinationNumberList = Measure.getCoordinationNumberList( corLabMap )
np.savetxt( (sample + 'cnList.txt'), coordinationNumberList, delimiter=',')
'''
sample = '2QR-1500N-2'
fileNameLoc = '/home/eg/codes/pacOutput/2QR-1500N-2/corLabMap.tiff'
print('Checking for ' + sample)
corLabMap = tf.imread( fileNameLoc ).astype( 'uint32' )
coordinationNumberList = Measure.getCoordinationNumberList( corLabMap )
np.savetxt( (sample + 'cnList.txt'), coordinationNumberList, delimiter=',')
'''
sample = 'OGF-0N'
fileNameLoc = '/home/eg/codes/pacOutput/OGF-0N/corLabMap.tiff'
print('Checking for ' + sample)
corLabMap = tf.imread( fileNameLoc ).astype( 'uint32' )
coordinationNumberList = Measure.getCoordinationNumberList( corLabMap )
np.savetxt( (sample + 'cnList.txt'), coordinationNumberList, delimiter=',')

sample = 'OGF-100N'
fileNameLoc = '/home/eg/codes/pacOutput/OGF-100N/corLabMap.tiff'
print('Checking for ' + sample)
corLabMap = tf.imread( fileNameLoc ).astype( 'uint32' )
coordinationNumberList = Measure.getCoordinationNumberList( corLabMap )
np.savetxt( (sample + 'cnList.txt'), coordinationNumberList, delimiter=',')

sample = 'OGF-500N'
fileNameLoc = '/home/eg/codes/pacOutput/OGF-500N/corLabMap.tiff'
print('Checking for ' + sample)
corLabMap = tf.imread( fileNameLoc ).astype( 'uint32' )
coordinationNumberList = Measure.getCoordinationNumberList( corLabMap )
np.savetxt( (sample + 'cnList.txt'), coordinationNumberList, delimiter=',')

sample = 'OGF-1500N-2'
fileNameLoc = '/home/eg/codes/pacOutput/OGF-1500N-2/corLabMap.tiff'
print('Checking for ' + sample)
corLabMap = tf.imread( fileNameLoc ).astype( 'uint32' )
coordinationNumberList = Measure.getCoordinationNumberList( corLabMap )
np.savetxt( (sample + 'cnList.txt'), coordinationNumberList, delimiter=',')


sample = 'OTC-0N'
fileNameLoc = '/home/eg/codes/pacOutput/OTC-0N/corLabMap.tiff'
print('Checking for ' + sample)
corLabMap = tf.imread( fileNameLoc ).astype( 'uint32' )
coordinationNumberList = Measure.getCoordinationNumberList( corLabMap )
np.savetxt( (sample + 'cnList.txt'), coordinationNumberList, delimiter=',')

sample = 'OTC-500N'
fileNameLoc = '/home/eg/codes/pacOutput/OTC-500N/corLabMap.tiff'
print('Checking for ' + sample)
corLabMap = tf.imread( fileNameLoc ).astype( 'uint32' )
coordinationNumberList = Measure.getCoordinationNumberList( corLabMap )
np.savetxt( (sample + 'cnList.txt'), coordinationNumberList, delimiter=',')

sample = 'OTC-1500N'
fileNameLoc = '/home/eg/codes/pacOutput/OTC-1500N/corLabMap.tiff'
print('Checking for ' + sample)
corLabMap = tf.imread( fileNameLoc ).astype( 'uint32' )
coordinationNumberList = Measure.getCoordinationNumberList( corLabMap )
np.savetxt( (sample + 'cnList.txt'), coordinationNumberList, delimiter=',')


sample = 'OTC-MD-7D50-0N'
fileNameLoc = '/home/eg/codes/pacOutput/OTC-MD-0N-7D50/corLabMap.tiff'
print('Checking for ' + sample)
corLabMap = tf.imread( fileNameLoc ).astype( 'uint32' )
coordinationNumberList = Measure.getCoordinationNumberList( corLabMap )
np.savetxt( (sample + 'cnList.txt'), coordinationNumberList, delimiter=',')

sample = 'OTC-MD-7D50-50N'
fileNameLoc = '/home/eg/codes/pacOutput/OTC-MD-50N-7D50/corLabMap.tiff'
print('Checking for ' + sample)
corLabMap = tf.imread( fileNameLoc ).astype( 'uint32' )
coordinationNumberList = Measure.getCoordinationNumberList( corLabMap )
np.savetxt( (sample + 'cnList.txt'), coordinationNumberList, delimiter=',')

sample = 'OTC-MD-7D50-500N'
fileNameLoc = '/home/eg/codes/pacOutput/OTC-MD-500N-7D50/corLabMap.tiff'
print('Checking for ' + sample)
corLabMap = tf.imread( fileNameLoc ).astype( 'uint32' )
coordinationNumberList = Measure.getCoordinationNumberList( corLabMap )
np.savetxt( (sample + 'cnList.txt'), coordinationNumberList, delimiter=',')

sample = 'OTC-MD-7D50-1500N'
fileNameLoc = '/home/eg/codes/pacOutput/OTC-MD-1500N-7D50/corLabMap.tiff'
print('Checking for ' + sample)
corLabMap = tf.imread( fileNameLoc ).astype( 'uint32' )
coordinationNumberList = Measure.getCoordinationNumberList( corLabMap )
np.savetxt( (sample + 'cnList.txt'), coordinationNumberList, delimiter=',')
'''
