# -*- coding: utf-8 -*-
"""
Outline:
This a part of the PAC code; this module is for visualization of the data and production of publication quality images

Features: 
    1. 2D cross sections of CT data
    2. 2D plots of X and Y parameters
    3. Rose diagrams of 2D orientation data
    4. Lambert projection plots of 3D orienation data
    5. 3D orientation data plots
    6. Intensity histogram plotter

References: 
    [1]

"""

# %% Importing libraries





#%%

class LemmeC:

    # Initialize
    def __init__(self):
        print("Visualizer activated")
    
    def plotGLI(self, arrayGLI):
        fig, (ax1, ax2, ax3) = plt.subplots(1,3) 
        ax1.imshow(arrayGLI[0], cmap='gray')
        ax2.imshow(arrayGLI[arrayGLI.shape[0]//2], cmap='gray')
        ax3.imshow(arrayGLI[arrayGLI.shape[0]-1], cmap='gray')          
    