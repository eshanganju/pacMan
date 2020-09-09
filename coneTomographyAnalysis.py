'''
Code for the analysis of 3D tomo data for samples collected around the cone penetrometer

written by @eg

Steps:
    Initialize the samples variables
    Read the files
    Extract subregions from the 3D tomo data - 6D50 edge length
    Compute particle sizes of the sand samples - user threshold after initial guess from OTSU
    Get particel aspect ration from PCA length ratios
    Get particle size distribution
    Calculate relative breakage parameter
    ---
    Get contact normals
    Get fabric tensors
    Get fabric projection plots

'''

import numpy as np
import skimage.external.tifffile as tf
import time

from pac import Reader
from pac import Segment
from pac import Measure
from pac import Plot

# Each test has three samples - top mid and bottom
# '2QR25', '2QR50', '2QR90','OGF25', 'OGF50', 'OGF90','OTC25', 'OTC50', 'OTC90'
scanData = ['2QR25']

for sample in scanData:

    # Initialize sample variables------------------------------------------------------------------------------------------------------*
    #   Input file location (3 for top, mid, and tip)
    #   Output file location
    #   Resolution
    #   Subregion size
    #   Shouler reference location in the top scan - this is a point which allows measurement of distances
    #   Dist between top mid and tip bottom left corner
    #   Starting subregion locations in each scan

    # 2QR sample-------------------*
    if sample == '2QR25':
        iflTop = '/home/eg/codes/pacInput/2QR_25_top/'
        oflTop = '/home/eg/codes/pacOutput/cone/2QR_25_top/'

        iflMid = '/home/eg/codes/pacInput/2QR_25_mid/'
        oflMid = '/home/eg/codes/pacOutput/cone/2QR_25_mid/'

        iflTip = '/home/eg/codes/pacInput/2QR_25_tip/'
        oflTip = '/home/eg/codes/pacOutput/cone/2QR_25_tip/'

        d50Ini = 0.73                                   # Initial D50 in mm

        resolutionTop = 0.01230                         # Resolution in mm/voxel
        resolutionMid = 0.01230                         # Resolution in mm/voxel
        resolutionTip = 0.01230                         # Resolution in mm/voxel

        subregionSize = 6*d50Ini                        # Subregion size is in mm (6 x D50)
        sholderRefLoc = np.array([513,317,392])         # Location of sholder in voxel units

        zTopMid = 0                                     # Z distance in mm between the lower left corners of Top and Mid scan
        yTopMid = 0                                     # Y distance in mm between the lower left corners of Top and Mid scan
        xTopMid = 0                                     # X distance in mm between the lower left corners of Top and Mid scan

        zTopTip = 0                                     # Z distance in mm between the lower left corners of Top and Tip scan
        yTopTip = 0                                     # Y distance in mm between the lower left corners of Top and Tip scan
        xTopTip = 0                                     # X distance in mm between the lower left corners of Top and Tip scan

        # Subregion top left corner in pixel units for the top scan location (ZZ, YY, XX)
        subregionLocTop = np.array( [ [],\
                                      [],\
                                      [],\
                                      [],\
                                      [],\
                                      [],\
                                      [] ] )

        # Subregion top left corner in pixel units for the mid scan location (ZZ, YY, XX)
        subregionLocMid = np.array( [ [],\
                                      [],\
                                      [],\
                                      [],\
                                      [],\
                                      [],\
                                      [] ] )

        # Subregion top left corner in pixel units for the tip scan location (ZZ, YY, XX)
        subregionLocTip = np.array( [ [],\
                                      [],\
                                      [],\
                                      [],\
                                      [],\
                                      [],\
                                      [] ] )

        # Single 3D matrix with details of the subregions
        ifl = np.array( [ iflTop , iflMid , iflTip ] )
        ofl = np.array( [ oflTop , oflMid , oflTip ] )
        resolution = np.array( [ resolutionTop, resolutionMid, resolutionTip ] )
        topMid = np.array( [ zTopMid , yTopMid , xTopMid ] )
        topTip = np.array( [ zTopTip , yTopTip , xTopTip ] )
        topDist = np.array( [ [ 0 , 0 , 0 ] , topMid , topTip ]  )
        subregionLoc = np.zeros( ( 3 , subregionLocTop.shape[0] , subregionLocTop.shape[1] ) )
        subregion[0] = subregionLocTop
        subregion[1] = subregionLocMid
        subregion[2] = subregionLocTip

    if sample == '2QR50':
        print('code goes plop...')

    if sample == '2QR90':
        print('code goes plop...')

    # OGF sample-------------------*
    if sample == 'OGF25':
        print('code goes plop...')

    if sample == 'OGF50':
        print('code goes plop...')

    if sample == 'OGF90':
        print('code goes plop...')

    # OTC sample-------------------*
    if sample == 'OTC25':
        print('code goes plop...')

    if sample == 'OTC50':
        print('code goes plop...')

    if sample == 'OTC90':
        print('code goes plop...')

    # Analysis--------------------------------------------------------------------------------------------------------------------------*
    # Analysis carried out for each section (top mid and top)
    for section in range(0,3):
        # Read the files----------------------------------------------------------------------------------------------------------------*
        #   Read one scan at a time and only that much in Z direction that is large enough to have the subregion size
        #   Extract subregions in a loop till all the subregions needed are analyzed. - no fabric so the gradation analysis should be faster

        # Get particle size distribution for subregions---------------------------------------------------------------------------------*
        #   Take subregions from the file and analyze for particle size distribution - Bin - EDM - WS - GSD
        #   No physical threhold - obtain intial guess form OTSU for threshold and then iterate till:
        #       1. Capture the largest particle size (using Feret diameter) same as the original gradation
        #       2. Visally no errors in segmentation - over or under segmentation
        #   Save labelled map and particle size distribution

        # Get from particle size distribution:------------------------------------------------------------------------------------------*
        #   Relative breakage parameter
        #   Morphology distribution
        #   Save relative breakage paramters and Morphology

        # Get from labelled map:--------------------------------------------------------------------------------------------------------*
        #   Contact normal distribution
        #   Fabric tensors
        #   Coordination number distribution
        #   Save fabric data




