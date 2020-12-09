'''
Code for the analysis of 3D tomo data for samples collected around the cone penetrometer

written by @eg

Steps:
    (01) Initialize the samples variables
    (02) Read the files and invert if needed
    (03) Extract subregions from the 3D tomo data - 6D50 edge length
    (04) Compute particle sizes of the sand samples - user threshold after initial guess from OTSU
    (05) Get particel aspect ration from PCA length ratios
    (06) Get particle size distribution
    (07) Calculate relative breakage parameter
    (08) Get contact normals
    (09) Get fabric tensors
    (10) Get fabric projection plots

'''

from pac import Reader
from pac import Segment
from pac import Measure
from pac import Plot

# Each test has three samples - top mid and bottom
scanData = ['2QR25','2QR25','2QR50','2QR90','OGF25','OGF50','OGF90','OTC25','OTC50','OTC90']
invertData = True
imageDataDump = ''

'''

'''
