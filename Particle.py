# -*- coding: utf-8 -*-

"""
Outline:

---
Features:


---
References: 


"""

# %% Importing libraries

import numpy as np



#%%

class Particle:

    def __init__(self, num, numVox, locData):
        self.index = num
        self.volume = numVox
        self.locationData = locData
        print("Particle No.",self.index,"created")