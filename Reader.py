# -*- coding: utf-8 -*-

"""
Purpose: 
    This class will read the data files from the local storage
    Convert the GLI from what is on the data to 0-1 range
"""

# General
import numpy as np
import skimage.external.tifffile as tiffy

# %% Classes
class Reader:

    def __init__(self):
        print("Reader activated")
    
    def imageRead(self, fileName):
        data = tiffy.imread(fileName)
        return data 



 
