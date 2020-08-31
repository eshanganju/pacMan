'''
Plotting
'''
import matplotlib.pyplot as plt
import spam
import numpy as np
import spam.plotting as splt
'''
Objectives are
    0. Cross-section plots
    1. Plot GSD
    2. Contact vector representations
        2.1 2D rose plots
        2.2 3D rose plots - blob plots
        2.3 Equal area projection plots
'''

VERBOSE = True

def centerCrossSection( labelledMap, colorMap = 'rainbow'):
    '''
    Show the first and last points of the axis
    Add a scale on the color map
    Make scale so that the image does not stretch
    '''

def grainSizeDistribution (gsd0, gsd1, gsd2=np.zeros((1,2)), gsd3=np.zeros((1,2)), gsd4=np.zeros((1,2)), gsd5=np.zeros((1,2)), gsd6=np.zeros((1,2)), xmax = 10, xmin = 0.001):
    '''
    plots semilogx grain size distribution
    add additional  gradations if passed
    add a legend with the gradations, ask user for legend labels
    '''
    plt.figure()
    plt.semilogx( gsd0[ :,0 ] , gsd0[ : , 1 ] , label='0' )
    plt.semilogx( gsd1[ :,0 ] , gsd1[ : , 1 ] , label='1' )
    if gsd2.sum() != 0: plt.semilogx( gsd2[ :,0 ] , gsd2[ : , 1 ] , label='2' )
    if gsd3.sum() != 0: plt.semilogx( gsd3[ :,0 ] , gsd3[ : , 1 ] , label='3' )
    if gsd4.sum() != 0: plt.semilogx( gsd4[ :,0 ] , gsd4[ : , 1 ] , label='4' )
    if gsd5.sum() != 0: plt.semilogx( gsd5[ :,0 ] , gsd5[ : , 1 ] , label='5' )
    if gsd6.sum() != 0: plt.semilogx( gsd6[ :,0 ] , gsd6[ : , 1 ] , label='6' )
    plt.xlim( 0.01 , 10 )
    plt.ylim( 0 , 100 )
    plt.legend()
    plt.show()

def rosePlot():
    '''
    '''

def blobPlot():
    '''
    '''

def equalAreaProjection( contactSummaryTable , makeYVertical = True, projectionOption='EA', plotToPresentOption='both', nr=9,xSize=8,ySize=4):
    print('\nPlotting the equal area projection of the contact normals')
    if makeYVertical == True: contactSummaryTable[:,[2,3]] = contactSummaryTable[:,[3,2]]

    orientations = contactSummaryTable[ : , 2 : 5 ]

    if projectionOption != 'EA': projectionUser = "equidistant"
    else: projectionUser = "lambert"

    if plotToPresentOption == 'points': plotToPresentUser = "points"
    elif plotToPresentOption == 'bin': plotToPresentUser = "bin"
    else : plotToPresentUser = "both"

    #maxBinVal = int(input('Enter max bin val: '))
    #minBinVal = int(input('Enter min bin val: '))

    maxBinVal = 10
    minBinVal = 0

    #normalize = input('Normalize bin([y]/n)?')
    normalize ='n'
    if normalize == 'n': normBinning = False
    else: normBinning = True

    print('Plotting results...')
    splt.plotOrientations( orientations_zyx = orientations,
                           projection = projectionUser,
                           plot = plotToPresentUser,
                           binValueMin = minBinVal,
                           binValueMax = maxBinVal,
                           binNormalisation = False,
                           numberOfRings = nr,
                           pointMarkerSize = 2,
                           cmap = plt.cm.gist_yarg,
                           title = "",
                           saveFigPath = None,
			   figx = xSize,
			   figy = ySize)
    plt.show()


