'''This code takes the feret psd and the value of particle volume loss and calculates the relative breakage value

The relative breakage is calculated using the Hardine (1985) and Einav (2007) defintions. 

Since in the analysis, we also remove small particles that are smaller than 1000 voxels in volume, we correct the gsd to take 
in to account the removed particle

In the case of Einav(2007), the ultimate grain size distribution is calculated assuming a fractal dimension of 2.6 and the 
dmax is assumed to be the smallest dmax for the inital, ultimate, and curent particle size distribution.

If the percentage passing for a particle size "i" in the uncorrected PSD is "p_uncorr_i" and the percentage of particles 
removed is "p_small," then the corrected percenage passing "p_corr_i" is computed as:
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

f = open( ouputLocation + 'breakageHardinDataSummary.csv',"w+" ) 
f.write('Name,PLoss,Np,BP,BC,BR')

# g = open( ouputLocation + 'breakageEinavDataSummary.csv',"w+" ) 
# g.write('Name,PLoss,Np,BP,BC,BR')

count = 1
countAll = len(fileList)

for psdFile in fileList:

	# Notify
	print('Analyzing ' + str(round(count)) + '/' + str(round(countAll)))
	
	# get psd
	psdComplete = np.loadtxt(psdFile, delimiter=',')
	psdOnly = psdComplete[:,-2:]

	# get percentage loss
	psdFileName = psdFile.split('/')[-1]
	dataName = psdFileName.split('-')
	plossFileName = pLossLocation + dataName[0] + '-' + dataName[1] + '-volumeLossSmallPtcl.txt'
	lossVal = float(np.genfromtxt(plossFileName,delimiter=',',skip_header=4,max_rows=1))

	# Update psd using percentage loss
	psdOnlyUpdate = np.zeros_like(psdOnly)
	psdOnlyUpdate[:,0] = psdOnly[:,0]
	psdOnlyUpdate[:,1] = ( psdOnly[:,1] + lossVal )/( 100 + lossVal) *100
	psdUpdateFileName = ouputLocation + dataName[0] + '-' + dataName[1] + '-gsdUpdate.csv'


	# Compute Hardine Br using updated pSD
	sandName = dataName[0].split('_')[0]
	
	if sandName == 'OTC': 
		psdOrig = np.loadtxt('/home/eg/codes/pacInput/originalGSD/otcOrig.csv',delimiter=',')
		psdOnlyUpdate[:,0] = psdOnlyUpdate[:,0] - 0.11 # offset - grey levels threshold issue?
	
	if sandName == 'OGF': 
		psdOrig = np.loadtxt('/home/eg/codes/pacInput/originalGSD/ogfOrig.csv',delimiter=',')
	
	if sandName == '2QR': 
		psdOrig = np.loadtxt('/home/eg/codes/pacInput/originalGSD/2qrOrig.csv',delimiter=',')

	# Save updated psd
	np.savetxt(psdUpdateFileName,psdOnlyUpdate,delimiter=',')

	# Save data
	numPtcl = psdOnlyUpdate.shape[0]
	
	# Hardin:
	if True:
		potBreakage, curBreakage, relBreakage = Measure.getRelativeBreakageHardin( psdOriginal=psdOrig, 
																				   psdCurrent=psdOnlyUpdate, 
																				   smallSizeLimit=0.075, 
																				   saveData=True, 
																				   sampleName= (dataName[0] + '-' + dataName[1]), 
																				   outputDir=ouputLocation )

		maxSize = psdOrig[ :,0 ].max()
		
		psdOnlyUpdate = np.append( psdOnlyUpdate, np.array( [maxSize , 100.0] ).reshape( 1, 2 ), 0 )

		minSize = 0.075
		psdClipped = psdOnlyUpdate[np.where(psdOnlyUpdate[:,0] > minSize)]
		minpp = min(psdClipped[:,1])
		val = np.array([minSize,minpp]).reshape(1,2)

		psdOnlyUpdate = np.concatenate((val,psdClipped), axis=0)


		plt.figure()
		plt.semilogx( psdOrig[ :,0 ] , psdOrig[ : , 1 ] , label='Original' )
		plt.semilogx( psdOnlyUpdate[ :,0 ] , psdOnlyUpdate[ : , 1 ] , label=(dataName[0] + '-' + dataName[1]) )
		plt.xlim( 0.01 , 10 )
		plt.ylim( 0 , 100 )
		plt.legend()
		plt.draw()
		plt.savefig( ouputLocation+ (dataName[0] + '-' + dataName[1]) + '-hrdPsd.tiff', transparent = False, dpi = 100 )
		plt.close()


		f.write( '\n' \
			     + str(dataName[0]+ '-' + dataName[1]) + ',' \
				 + str(lossVal) + ',' \
				 + str(numPtcl) + ',' \
				 + str(potBreakage) + ',' \
				 + str(curBreakage) + ',' \
				 + str(relBreakage) )

	# Einav
	if False:
		smallSizeLimitForEinav=0.0001
		uPSD, potBreakage, curBreakage, relBreakage = Measure.getRelativeBreakageEinav( psdOriginal=psdOrig, 
																						psdCurrent=psdOnlyUpdate, 
																						fracDim=2.6, 
																						smallSizeLimit=smallSizeLimitForEinav,
																						saveData=True, 
																						sampleName= (dataName[0] + '-' + dataName[1]), 
																						outputDir=ouputLocation )

		maxSize = psdOrig[ :,0 ].max()
		
		psdOnlyUpdate = np.append( psdOnlyUpdate, np.array( [maxSize , 100.0] ).reshape( 1, 2 ), 0 )

		minSize = smallSizeLimitForEinav
		psdClipped = psdOnlyUpdate[np.where(psdOnlyUpdate[:,0] > minSize)]
		minpp = 0
		val = np.array([minSize,minpp]).reshape(1,2)
		psdOnlyUpdate = np.concatenate((val,psdClipped), axis=0)

		plt.figure()
		plt.semilogx( psdOrig[ : , 0 ] , psdOrig[ : , 1 ] , label='Original' )
		plt.semilogx( psdOnlyUpdate[ : , 0 ] , psdOnlyUpdate[ : , 1 ] , label=(dataName[0] + '-' + dataName[1]) )
		plt.semilogx( uPSD[: , 0],uPSD[ : , 1 ], label='Ultimate')
		plt.xlim( minSize , 10 )
		plt.ylim( 0 , 100 )
		plt.legend()
		plt.draw()
		plt.savefig( ouputLocation+ (dataName[0] + '-' + dataName[1]) + '-envPsd.tiff', transparent = False, dpi = 100 )
		plt.close()


		g.write( '\n' \
			     + str(dataName[0]+ '-' + dataName[1]) + ',' \
				 + str(lossVal) + ',' \
				 + str(numPtcl) + ',' \
				 + str(potBreakage) + ',' \
				 + str(curBreakage) + ',' \
				 + str(relBreakage) )

	count = count + 1


f.close()
# g.close()