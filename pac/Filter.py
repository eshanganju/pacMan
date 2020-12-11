'''
Filter class
'''

import tifffile as tiffy
from skimage import restoration
import matplotlib.pyplot as plt
import numpy as np
import time

def filterDenoiseNlm(gli, gliMax, outputDir, sandName='sampleX'):
    """
    Description
        This function takes in a unfiltered grey level intensity map and loops over parameters
        till it outputs a filtered map that is filtered using the non local means filter.
        The non local means filter removes noise with minimal effect on grey-level edges.
    Parameters
        gli (np array): the ct scan data in the form of a NP array.
        gliMax (int): the maximum integer value in the plot. (can be set to the max of the image or max of the dtype).
        outputDir (str): path to storage of the filtered image.
        sandName (str): name of the sand used for the analysis, default to sample X.
    Returns
        filtered image in the form of a numpy array.
    """
    inputImage = gli[ gli.shape[ 0 ] // 2 ]
    filteredImage = np.zeros_like( inputImage )
    noiseRemoved = np.zeros_like( filteredImage )

    # Local variables for the entire map
    inputMap = gli
    filteredMap = np.zeros_like( inputMap )

    # Histograms
    numHistPts = inputImage.shape[ 0 ] * inputImage.shape[ 1 ]
    numHistPtsMap = inputMap.shape[ 0 ] * inputMap.shape[ 1 ] * inputMap.shape[ 2 ]
    his_x = np.arange( 0, gliMax + 1, 1 )


    # The SD of the input noisy image
    sigmaEstimateImage = np.mean(restoration.estimate_sigma(inputImage, multichannel=True))
    sigmaEstimateMap = np.mean(restoration.estimate_sigma(inputMap, multichannel=True))


    # Inital settings chosen based on parameteric study
    patchSize = 3           # Patch size: 
    patchDistance = 5       # Patch distance: 
    hVal = 450              # h-value: 

    print('\n\nIntial parameters for filter----------------------*\n')
    print('Patch size: ', patchSize)
    print('Patch distance: ', patchDistance)
    print('Cut-off pixel intensity: ', round(hVal))


    # Loop for checking filtering parameters; run on CS of data
    runSliceFilter = True
    while runSliceFilter == True:
        print('\n\nNew filtering loop started on central cross-section--------*\n')
        start_time = time.time()
        filteredImage = restoration.denoise_nl_means(inputImage, patch_size = patchSize, patch_distance = patchDistance, h = hVal, multichannel = False, fast_mode = True, sigma = sigmaEstimateImage)
        noiseRemoved = abs( inputImage - filteredImage )
        listNoisy = inputImage.reshape( ( numHistPts, 1 ) )
        listClean = filteredImage.reshape( ( numHistPts, 1 ) )

        # Histogram
        histNoisy = np.histogram(listNoisy, bins = gliMax + 1, range = ( 0, gliMax ) )
        histClean = np.histogram(listClean, bins = gliMax + 1, range = ( 0, gliMax ) )
        histNoisy = histNoisy / ( histNoisy[0].sum() ) * 100
        histClean = histClean / ( histClean[0].sum() ) * 100

        # Figures
        plt.figure()
        plt.imshow(noiseRemoved, cmap= 'Greys_r' )
        plt.draw() # draw the plot
        nameofTempFile = outputDir + 'NoiseRemoved-CenterSlice.png'
        plt.savefig(nameofTempFile)
        plt.pause( 5 ) # show it for 5 seconds
        plt.close()

        # Figures
        plt.figure()
        plt.imshow(filteredImage, cmap= 'Greys_r' )
        plt.draw() # draw the plot
        nameofTempFile = outputDir + 'filteredImage-CenterSlice.png'
        plt.savefig(nameofTempFile)
        plt.pause( 5 ) # show it for 5 seconds
        plt.close()

        # Figures - Histogram
        plt.figure()
        plt.plot(his_x, histNoisy[ 0 ])
        plt.plot(his_x, histClean[ 0 ])
        plt.grid()
        plt.draw() # draw the plot
        nameofTempFile = outputDir + 'histogram-CenterSlice.png'
        plt.savefig(nameofTempFile)
        plt.pause( 5 ) # show it for 5 seconds
        plt.close()

        timeTakenThisLoop = ( time.time() - start_time )
        print( "\n--- Time taken: %s seconds ---" %round( timeTakenThisLoop ) )

        answer = input( "\n\nCheck files - are filter parameters suitable ([y]/n)?:" )
        if answer == 'n':
            print( "Enter new parameters: \n\n" )
            patchSize = int( input( "Patch size (pixel): " ) )
            patchDistance = int( input( "Patch distance (pixel):" ) )
            hVal = float( input( "Cut-off pixel intensity: " ) )
        else:
            runSliceFilter = False


    # Filtering entire grey level map with chosen parameters
    print( '\n\nFilter for entire grey level map started--------*' )
    print( 'This take a lot of time...' )
    start_time = time.time()
    filteredMap = restoration.denoise_nl_means(inputMap, patch_size=patchSize, patch_distance=patchDistance, h=hVal, multichannel=False, fast_mode=True, sigma=sigmaEstimateMap)
    listNoisy = inputMap.reshape( ( numHistPtsMap, 1 ) )
    listClean = filteredMap.reshape( ( numHistPtsMap, 1 ) )

    # Histogram
    histNoisy = np.histogram( listNoisy, bins = gliMax + 1, range = ( 0, gliMax ) )
    histClean = np.histogram( listClean, bins = gliMax + 1, range = ( 0, gliMax ) )
    histNoisy = histNoisy / ( histNoisy[ 0 ].sum() ) * 100
    histClean = histClean / ( histClean[ 0 ].sum() ) * 100

    # Figures - Histogram
    plt.figure()
    plt.plot(his_x, histNoisy[ 0 ])
    plt.plot(his_x, histClean[ 0 ])
    plt.grid()
    plt.draw() # draw the plot
    volumeHistogramName = outputDir + sandName+'GLIhistogram-volume.png'
    plt.savefig(volumeHistogramName)
    plt.pause( 5 ) # show it for 5 seconds
    plt.close()

    print( '.\n.\n.\nFilter for entire grey level map completed--------*' )
    timeTakenThisLoop = ( time.time() - start_time )
    print( "\n--- Time taken: %s minutes ---\n" % round( timeTakenThisLoop // 60 ) )

    # Updating filtered parameters
    filterParameterFileName = outputDir + sandName + '-NLMParameters.txt'
    f = open(filterParameterFileName,"w+")
    f.write("NLM filter parameters---------------------*\n")
    f.write("Patch size = %f\n" % patchSize)
    f.write("Patch distance = %f\n" % patchDistance)
    f.write("Cut-off intensity = %f\n" % hVal)
    f.write("Estimated sigma - CS = %f\n" % sigmaEstimateImage)
    f.write("Estimated sigma - Map = %f\n" % sigmaEstimateMap)
    f.write("Time taken for filter (s) = %f\n" % round(timeTakenThisLoop))
    f.close()

    return filteredMap


