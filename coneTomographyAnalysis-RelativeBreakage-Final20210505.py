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

inputFileLoc = '/home/eg/codes/pacInput/gradationsConeCrushFinal20210505/psd/'
pLossLocation = '/home/eg/codes/pacInput/gradationsConeCrushFinal20210505/volLoss/'

ouputLocation = '/home/eg/codes/pacOutput/cone/finalFolder/gsdAnalysis/'

searchString = inputFileLoc + '*csv'
fileList = glob.glob(searchString)

f = open( ouputLocation + 'breakageDataSummary.csv',"w+" ) 
f.write('Name,PLoss,Np,BP,BC,BR')

for gsdFile in fileList:
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

	# Save updated gsd
	np.savetxt(gsdUpdateFileName,gsdOnlyUpdate,delimiter=',')

	# Compute Hardine Br using updated GSD
	sandName = dataName[0].split('_')[0]
	if sandName == 'OTC': psdOrig = np.loadtxt('/home/eg/codes/pacInput/originalGSD/otcOrig.csv',delimiter=',')
	if sandName == 'OGF': psdOrig = np.loadtxt('/home/eg/codes/pacInput/originalGSD/ogfOrig.csv',delimiter=',')
	if sandName == '2QR': psdOrig = np.loadtxt('/home/eg/codes/pacInput/originalGSD/2qrOrig.csv',delimiter=',')
	
	potBreakage, curBreakage, relBreakage = Measure.getRelativeBreakageHardin( psdOriginal=psdOrig, 
																			   psdCurrent=gsdOnlyUpdate, 
																			   smallSizeLimit=0.075, 
																			   saveData=True, 
																			   sampleName= (dataName[0] + '-' + dataName[1]), 
																			   outputDir=ouputLocation )

	plt.figure()
	plt.semilogx( psdOrig[ :,0 ] , psdOrig[ : , 1 ] , label='Original' )
	plt.semilogx( gsdOnlyUpdate[ :,0 ] , gsdOnlyUpdate[ : , 1 ] , label=(dataName[0] + '-' + dataName[1]) )
	plt.xlim( 0.01 , 10 )
	plt.ylim( 0 , 100 )
	plt.legend()
	plt.draw()
	plt.savefig( ouputLocation+ (dataName[0] + '-' + dataName[1]) + '-gsd.tiff', transparent = False, dpi = 100 )
	plt.close()


	# Save data
	numPtcl = gsdOnlyUpdate.shape[0]
	
	f.write( '\n' \
		     + str(dataName[0]) + ',' \
			 + str(lossVal) + ',' \
			 + str(numPtcl) + ',' \
			 + str(potBreakage) + ',' \
			 + str(curBreakage) + ',' \
			 + str(relBreakage) )

f.close()