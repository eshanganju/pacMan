# -*- coding: utf-8 -*-

"""
Created on  2019-10-15
Level set implementation in pyhton
"""

import numpy as np
import matplotlib.pyplot as plt
import skimage.external.tifffille as tiffy

fglm = tiffy.imread('FilteredBox4B.tiff')
bi = tiffy.imread('Bin.tiff')
mrk = tiffy.imread('Markers.tiff')
seg = tiffy.imread('watershedSegmentation.tiff')



