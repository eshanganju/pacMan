'''
Code for the analysis of 3D tomo data for samples collected around the cone penetrometer

written by @eg
'''

import numpy as np
import skimage.external.tifffile as tf
import time

from pac import Reader
from pac import Segment
from pac import Measure
from pac import Plot

scanData = ['2QR25', '2QR50', '2QR90','OGF25', 'OGF50', 'OGF90','OTC25', 'OTC50', 'OTC90']

for sample in scanData:
    '''
    Initialize the samples variables
    Read the files
    `
    '''

