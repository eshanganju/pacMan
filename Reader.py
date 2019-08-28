# -*- coding: utf-8 -*-

"""
This class will read the data files from the local storage
Maybe even generate new files using Kalisphera - not sure yet
"""

# General
import numpy as np
import skimage.external.tifffile as tiffy


class Reader:

    def __init__(self):
        print("Reader activated")
    
    def imageRead(self, fileName):
        data = tiffy.imread(fileName)
        return data 



 
