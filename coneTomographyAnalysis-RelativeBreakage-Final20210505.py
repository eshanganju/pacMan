'''This code takes the feret psd and the value of particle volume loss and calculates the relative breakage value

The relative breakage is calculated using the Hardine (1985) definition. The reason for this is that the the largest particle in the subrgegion can be smaller than the largest particle in the original sample

Since in the analysis, we also remove small particles that are smaller than 1000 voxels in volume, we correct the gsd to take in to account the removed particle

If the percentage passing for a particle size "i" in the uncorrected PSD is "p_uncorr_i" and the percentage of particles removed is "p_small," then the corrected percenage passing "p_corr_i" is computed as:
	p_corr_i = ( p_uncorr_i + p_small )/( 100 + p_small ) * 100

'''

from pac import Measure
from pac import Read
import glob
import numpy as np
import matplotlib.pyplot as plt

# Read file list from the folder

# Loop over file list
	# Read percentage loss for corresponding file
	# Read uncorrected particle size distribution
	# Updated uncorrected particle size distribution

inputFileLoc = '/home/eg/codes/pacInput/gradationsConeCrushFinal20210505-2/psd/'
pLossLocation = '/home/eg/codes/pacInput/gradationsConeCrushFinal20210505-2/volLoss/'

ouputLocation = '/home/eg/codes/pacOutput/cone/finalFolder/gsdAnalysis/'

searchString = inputFileLoc + '*csv'
fileList = glob.glob(searchString)

f = open( ouputLocation + 'breakageDataSummary.csv',"w+" ) 
f.write('Name,PLoss,Np,BP,BC,BR')

count = 1
countAll = len(fileList)

for gsdFile in fileList:

	# Notify
	print('Analyzing ' + str(round(count)) + '/' + str(round(countAll)))
	
	# get gsd
	gsdComplete = np.loadtxt(gsdFile, delimiter=',')
	gsdOnly = gsdComplete[:,-2:]

	# get percentage loss
	gsdFileName = gsdFile.split('/')[-1]
	dataName = gsdFileName.split('-')
	plossFileName = pLossLocation + dataName[0] + '-' + dataName[1] + '-volumeLossSmallPtcl.txt'
	lossVal = float(np.genfromtxt(plossFileName,delimiter=',',skip_header=4,max_rows=1))

	# Update gsd using percentage loss
	gsdOnlyUpdate = np.zeros_like(gsdOnly)
	gsdOnlyUpdate[:,0] = gsdOnly[:,0]
	gsdOnlyUpdate[:,1] = ( gsdOnly[:,1] + lossVal )/( 100 + lossVal) *100
	gsdUpdateFileName = ouputLocation + dataName[0] + '-' + dataName[1] + '-gsdUpdate.csv'


	# Compute Hardine Br using updated GSD
	sandName = dataName[0].split('_')[0]
	
	if sandName == 'OTC': 
		psdOrig = np.loadtxt('/home/eg/codes/pacInput/originalGSD/otcOrig.csv',delimiter=',')
		gsdOnlyUpdate[:,0] = gsdOnlyUpdate[:,0] - 0.11
	
	if sandName == 'OGF': 
		psdOrig = np.loadtxt('/home/eg/codes/pacInput/originalGSD/ogfOrig.csv',delimiter=',')
	
	if sandName == '2QR': 
		psdOrig = np.loadtxt('/home/eg/codes/pacInput/originalGSD/2qrOrig.csv',delimiter=',')

	# Save updated gsd
	np.savetxt(gsdUpdateFileName,gsdOnlyUpdate,delimiter=',')

	# Save data
	numPtcl = gsdOnlyUpdate.shape[0]
	
	potBreakage, curBreakage, relBreakage = Measure.getRelativeBreakageHardin( psdOriginal=psdOrig, 
																			   psdCurrent=gsdOnlyUpdate, 
																			   smallSizeLimit=0.075, 
																			   saveData=True, 
																			   sampleName= (dataName[0] + '-' + dataName[1]), 
																			   outputDir=ouputLocation )

	maxSize = psdOrig[ :,0 ].max()
	
	gsdOnlyUpdate = np.append( gsdOnlyUpdate, np.array( [maxSize , 100.0] ).reshape( 1, 2 ), 0 )

	minSize = 0.075
	psdClipped = gsdOnlyUpdate[np.where(gsdOnlyUpdate[:,0] > minSize)]
	minpp = min(psdClipped[:,1])
	val = np.array([minSize,minpp]).reshape(1,2)

	gsdOnlyUpdate = np.concatenate((val,psdClipped), axis=0)


	plt.figure()
	plt.semilogx( psdOrig[ :,0 ] , psdOrig[ : , 1 ] , label='Original' )
	plt.semilogx( gsdOnlyUpdate[ :,0 ] , gsdOnlyUpdate[ : , 1 ] , label=(dataName[0] + '-' + dataName[1]) )
	plt.xlim( 0.01 , 10 )
	plt.ylim( 0 , 100 )
	plt.legend()
	plt.draw()
	plt.savefig( ouputLocation+ (dataName[0] + '-' + dataName[1]) + '-gsd.tiff', transparent = False, dpi = 100 )
	plt.close()


	f.write( '\n' \
		     + str(dataName[0]+ '-' + dataName[1]) + ',' \
			 + str(lossVal) + ',' \
			 + str(numPtcl) + ',' \
			 + str(potBreakage) + ',' \
			 + str(curBreakage) + ',' \
			 + str(relBreakage) )

	count = count + 1


f.close()