'''
Plotting
'''
import matplotlib.pyplot as plt
import spam

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

def grainSizeDistribution (gsd1, gsd2=None, gsd3=None, gsd4=None, gsd5=None, gsd6=None, xmax = 10, xmin = 0.001):
    '''
    plots semilogx grain size distribution
    add additional  gradations if passed
    add a legend with the gradations, ask user for legend labels
    '''
    plt.figure()
    plt.semilogx( gsd1[ :,0 ] , gsd1[ : , 1 ] , label='1' )
    if gsd2 != None : plt.semilogx( gsd2[ :,0 ] , gsd2[ : , 1 ] , label='2' )
    if gsd3 != None : plt.semilogx( gsd3[ :,0 ] , gsd3[ : , 1 ] , label='3' )
    if gsd4 != None : plt.semilogx( gsd4[ :,0 ] , gsd4[ : , 1 ] , label='4' )
    plt.xlim( 0.001 , 10 )
    plt.ylim( 0 , 100 )
    plt.legend()
    plt.show()

def rosePlot2D():
    '''
    '''

def rosePlot3D():
    '''
    '''

def equalAreaProjection( contactTable ):
    '''
    '''

