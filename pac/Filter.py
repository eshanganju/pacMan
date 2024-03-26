'''Filter module: runs a loop to clean out the image.
'''

import tifffile as tiffy
from skimage import restoration
import matplotlib.pyplot as plt
import numpy as np
import time

# This is to plot all the text when running the functions
VERBOSE = True
TESTING = True

def filterUsingNlm(gli, bitDepth=16, pSize=3, pDistance=7, hVal=None, saveImg=False, outputDir='', sampleName='sampleX', loop=True):
	"""This function takes in a unfiltered grey level and filters it using non-local means (nlm) filter
	This function loops over nlm parameters till user accepts the filtered image.

	The nlm filter removes noise with minimal effect on grey-level edges.

	Parameters
	----------
	gli : unsigned integer ndarray
		The XCT data in the form of a NP array

	bitDepth : int
		The bit depth of the XCT data

	pSize : unsigned integer 
		patch size for application of nlm filter

	pDist : integer
		Edge of volume over which the image is seacherd for mean

	hVal : integer
		cut off value for filter GLI

	saveImg : bool
		Should the fileterd image be saved?

	outputDir : string
		Location of output directory

	sampleName : string 
		Name of the sample used for the analysis, default to sampleX

	Return
	------
	filteredMap : unsigned integer ndarray
		Filtered image in the form of a numpy array.
	"""

	print('Starting non-local means filter')
	print('-------------------------------\n\n')

	gliMax = 2**bitDepth-1

	inputImage = gli[ gli.shape[ 0 ] // 2 ]
	filteredImage = np.zeros_like( inputImage )
	noiseRemoved = np.zeros_like( filteredImage )

	inputMap = gli
	filteredMap = np.zeros_like( inputMap )

	numHistPts = inputImage.shape[ 0 ] * inputImage.shape[ 1 ]
	numHistPtsMap = inputMap.shape[ 0 ] * inputMap.shape[ 1 ] * inputMap.shape[ 2 ]
	his_x = np.arange( 0, gliMax + 1, 1 )

	sigmaEstimateImage = restoration.estimate_sigma(inputImage)
	sigmaEstimateMap = restoration.estimate_sigma(inputMap)

	if hVal==None: hVal = 1.5*sigmaEstimateImage 

	print('Gaussian noise slice: ' + str(sigmaEstimateImage) + '\n')
	print('Gaussian noise Map: ' + str(sigmaEstimateMap) + '\n')

	noiseParameterFileName = outputDir+sampleName+'-gaussianNoise.txt'
	f = open(noiseParameterFileName,"w+")
	f.write("Noise-------*\n")
	f.write("Gaussian noise slice: %f\n" %sigmaEstimateImage)
	f.write("Gaussian noise map: %f\n" %sigmaEstimateMap)
	f.close()

	print('\n\nIntial parameters for filter----------------------*\n')
	print('Patch size: ', pSize)
	print('Patch distance: ', pDistance)
	print('Cut-off pixel intensity: ', round(hVal))
	print('--------------------------------------------------*\n')

	runSliceFilter = True
	while runSliceFilter == True:

		print('\n\nNew filtering loop started on central cross-section--------*\n')
		start_time = time.time()

		print('\n\nCurrent parameters for filter----------------------*\n')
		print('Patch size: ', pSize)
		print('Patch distance: ', pDistance)
		print('Cut-off pixel intensity: ', round(hVal))
		print('--------------------------------------------------*\n')

		filteredImage = restoration.denoise_nl_means( inputImage,
														patch_size=pSize,
														patch_distance=pDistance,
														preserve_range=True,
														h=hVal,
														fast_mode=True,
														sigma=sigmaEstimateMap )

		noiseRemoved = abs( inputImage - filteredImage )

		listNoisy = inputImage.reshape( ( numHistPts, 1 ) )
		listClean = filteredImage.reshape( ( numHistPts, 1 ) )
		
		"""
		histNoisy = np.histogram(listNoisy, bins = gliMax + 1, range = ( 0, gliMax ) )
		histClean = np.histogram(listClean, bins = gliMax + 1, range = ( 0, gliMax ) )
		histNoisy = histNoisy / ( histNoisy[0].sum() ) * 100
		histClean = histClean / ( histClean[0].sum() ) * 100
		np.savetxt(outputDir+sampleName+'-histSliceNoisy.csv',histNoisy[0],delimiter=',',header='%')
		np.savetxt(outputDir+sampleName+'-histSliceClean.csv',histClean[0],delimiter=',',header='%')
		"""

		plt.figure()
		plt.imshow(inputImage, cmap= 'Greys_r' )
		plt.draw() # draw the plot
		nameofTempFile = outputDir+sampleName+'-noisyImage-CenterSlice.tif'
		tiffy.imsave(nameofTempFile,inputImage)
		plt.pause( 1 ) # show it for 1 seconds
		plt.close()

		filterDetails = '-ps' + str(round(pSize)) + '-pd' + str(round(pDistance)) + '-h' + str(round(hVal))

		plt.figure()
		plt.imshow(filteredImage, cmap= 'Greys_r' )
		plt.draw() # draw the plot
		nameofTempFile = outputDir+sampleName+'-filteredImage-CenterSlice' + filterDetails + '.tif'
		tiffy.imsave(nameofTempFile,filteredImage)
		plt.pause( 1 ) # show it for 5 seconds
		plt.close()

		plt.figure()
		plt.imshow(noiseRemoved, cmap= 'Greys_r' )
		plt.draw() # draw the plot
		nameofTempFile = outputDir+sampleName+'-noiseRemoved-CenterSlice' + filterDetails + '.tif'
		tiffy.imsave(nameofTempFile,noiseRemoved)
		plt.pause( 1 ) # show it for 5 seconds
		plt.close()

		"""
		plt.figure()
		plt.plot(his_x, histNoisy[ 0 ])
		plt.plot(his_x, histClean[ 0 ])
		plt.grid()
		plt.draw() # draw the plot
		nameofTempFile = outputDir+sampleName+'-histogram-CenterSlice' + filterDetails + '.tif'
		plt.savefig(nameofTempFile)
		plt.pause( 1 ) # show it for 5 seconds
		plt.close()
		"""
		
		timeTakenThisLoop = ( time.time() - start_time )
		print( "\n--- Time taken: %s seconds ---" %round( timeTakenThisLoop ) )

		if loop == True:
			answer = input( "\n\nCheck files - are filter parameters suitable ([y]/n)?:" )
			if answer == 'n':
				print( "Enter new parameters: \n\n" )
				pSize = int( input( "Patch size (pixel): " ) )
				pDistance = int( input( "Patch distance (pixel):" ) )
				hVal = float( input( "Cut-off pixel intensity: " ) )
			else:
				runSliceFilter = False

		if loop == False: runSliceFilter = False


	print( '\n\nFilter for entire grey level map started--------*' )
	print( 'This take a lot of time...' )
	start_time = time.time()

	filteredMap = restoration.denoise_nl_means(inputMap,
												patch_size=pSize,
												patch_distance=pDistance,
												h=hVal,
												fast_mode=True,
												sigma=sigmaEstimateMap)

	listNoisy = inputMap.reshape( ( numHistPtsMap, 1 ) )
	listClean = filteredMap.reshape( ( numHistPtsMap, 1 ) )

	histNoisy = np.histogram( listNoisy, bins = gliMax + 1, range = ( 0, gliMax ) )
	histClean = np.histogram( listClean, bins = gliMax + 1, range = ( 0, gliMax ) )
	histNoisy = histNoisy / ( histNoisy[ 0 ].sum() ) * 100
	histClean = histClean / ( histClean[ 0 ].sum() ) * 100
	np.savetxt(outputDir+sampleName+'-histNoisy.csv',histNoisy[0],delimiter=',',header='%')
	np.savetxt(outputDir+sampleName+'-histClean.csv',histClean[0],delimiter=',',header='%')

	plt.figure()

	plt.plot(his_x,histNoisy[ 0 ])
	plt.plot(his_x,histClean[ 0 ])
	plt.grid()
	plt.draw() # draw the plot
	volumeHistogramName = outputDir+sampleName+'-GLIhistogram-volume.png'
	plt.savefig(volumeHistogramName)
	plt.pause( 5 ) # show it for 5 seconds
	plt.close()

	timeTakenThisLoop = ( time.time() - start_time )

	if VERBOSE: 
		print( '.\n.\n.\nFilter for entire grey level map completed--------*' )
		print( "\n--- Time taken: %s minutes ---\n" % round( timeTakenThisLoop // 60 ) )

	filterParameterFileName = outputDir+sampleName+'-NLMParameters.txt'
	f = open(filterParameterFileName,"w+")
	f.write( "NLM filter parameters---------------------*\n")
	f.write( "Patch size = %f\n" % pSize )
	f.write( "Patch distance = %f\n" % pDistance )
	f.write( "Cut-off intensity = %f\n" % hVal)
	f.write( "Estimated sigma - CS = %f\n" % sigmaEstimateImage )
	f.write( "Estimated sigma - Map = %f\n" % sigmaEstimateMap )
	f.write( "Time taken for filter (s) = %f\n" % round(timeTakenThisLoop ) )
	f.close()

	if saveImg == True: tiffy.imsave(outputDir+sampleName+'-filteredGLIMap.tif',filteredMap.astype('uint16'))

	return filteredMap

def _filterUsingMedian(gli, bitDepth=16, pSize=3, outputDir='',sampleName='',saveData=True,saveFilterParams=True,loop=True, slice='central'):
	"""


	Parameters
	----------

		slice: string; either 'central', 'upper', 'lower', 'allThree' to carry out intial filteration on a slice in the center, upper quartile, lower quartile, or all three positions.

	Return
	------
		filteredMap: ndArray; numpy array the same size as the input image

	"""
	filteredMap = gli

	# Print intial filtration settings

	# Start loop
		# Extract and filter slices of the image
		# Show the images and noise removed + noise removal parameters

		# Prepare image for filtration
	filrationIsGood = False

		# Loop if filtrationIsGood == False

		# When filtrationIsGood == True
			# Filter entire dataset in 3D using filter parameters
			# Save filter parameters if saveFilterParams == True

	# Return filtered dataset
	return filteredMap












