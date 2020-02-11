#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Purpose: 
    This class will carry out misc activities to make sure all function alright

Jeeves: 
    Keeps the location of input and output files
    
"""

# General
import numpy as np
import skimage.external.tifffile as tiffy

# %% Class
class Jeeves:

    def __init__(self):
        print( 'Activating Jeeves...' )
        self.sandName = ''
        self.inputFilesLocation = ''
        self.outputFilesLocation = ''
        self.superCubeEdgeLengthforREVAnalysis = float( 0 )
        self.superCubeCenterSlice = int( 0 )
        self.D50 = float( 0 )
        self.superCubeCenterRow = int( 0 ) # Y
        self.superCubeCenterCol = int( 0 ) # X
        self.tiffFileLocation = ''
        self.superCubeCalib = float( 0 )
        self.voidRatio = float(0)
        print( 'Jeeves activated' )

             
    def checkFolderLocations(self):
        
        defaultInputFilesLocation = 'C:/Users/eganj/gitHub/pacInput/'
        defaultOutputFilesLocation = 'C:/Users/eganj/gitHub/pacOutput/'       
        
        # Check default input: 
        inputLocationIsCorrect = input( 'By default, input files are supposed to be at: ' + defaultInputFilesLocation + '\nIs this coorect? [y/n]: ' )        
        if inputLocationIsCorrect.lower() == 'y':
            self.inputFilesLocation = defaultInputFilesLocation       
        else:
            self.inputFilesLocation = input( 'Enter new location of input files: ' )           

        # Check default output: 
        outputLocationIsCorrect = input( '\nBy default, output files are supposed to be at: ' + defaultOutputFilesLocation + '\nIs this coorect? [y/n]: ' )                 
        if outputLocationIsCorrect.lower() == 'y':
            self.outputFilesLocation = defaultOutputFilesLocation      
        else:
            self.outputFilesLocation = input('Enter new location of output files: ')

        print('\nFolder locations established...')
    
    
    def checkREVSampleDetails(self):
        
        # supercube calibration
        print( 'superCube default calibration is 0.01193 mm/px...' )
        changeSuperCubeCalib = input( 'Do you want to change this [y/n]?: ' )
        if changeSuperCubeCalib.lower() == 'y':
            self.superCubeCalib = float( input( 'Enter calibration (mm/px): ' ) )        
        else:
            self.superCubeCalib = 0.01193
                
        # Supercube size details: 
        print( 'superCube default edge length is 5.5 mm' )
        changeSuperCubeEdgeLength = input( 'Do you want to change this [y/n]?: ' )      
        if changeSuperCubeEdgeLength.lower() == 'y':
            self.superCubeEdgeLengthforREVAnalysis = float( input( 'Enter superCube edge length (mm): ' ) )
        else:
            self.superCubeEdgeLengthforREVAnalysis = 5.5
                
        # Sand sample details: 
        allSampleDetailsObtained = False       
        while allSampleDetailsObtained == False:
            sand = input( 'Enter sand name (ogf, 2qr, otc or new): ' )
            if sand.lower() == 'ofg':
                self.sandName = sand
                self.superCubeCenterSlice = round( ( 21+1005 ) / 2 )
                self.D50 = 0.62
                self.superCubeCenterRow = 450 # Y
                self.superCubeCenterCol = 506 # X
                self.voidRatio = 0.6348
                self.tiffFileLocation = self.inputFilesLocation+'OGF-0N'
                allSampleDetailsObtained = True    
            elif sand.lower() == 'otc':
                self.sandName = sand
                self.superCubeCenterSlice = round( ( 20+1006 ) / 2 )
                self.D50 = 0.72
                self.superCubeCenterRow = 450 # Y
                self.superCubeCenterCol = 506 # X
                self.voidRatio = 0.5404
                self.tiffFileLocation = self.inputFilesLocation+'OTC-0N'
                allSampleDetailsObtained = True                
            elif sand.lower() == '2qr':
                self.sandName = sand
                self.superCubeCenterSlice = round( ( 21+1005 ) / 2 )
                self.calib = 11930 / 1000 / 1000 # mm/px
                self.D50 = 0.71
                self.superCubeCenterRow = 450 # Y
                self.superCubeCenterCol = 506 # X
                self.voidRatio = 0.7337
                self.tiffFileLocation = self.inputFilesLocation+'2QR-0N'
                allSampleDetailsObtained = True
            elif sand.lower() == 'new':
                self.sandName = input('Enter sand name: ')
                self.superCubeCenterSlice = int( input( 'Enter the center slice number: ' ))
                self.calib = float( input( 'Enter the calibration (mm/px): ' ) )
                self.D50 = float( input( 'Enter the D50 of the soil (mm): ' ) )
                self.superCubeCenterRow = int( input( 'Enter the center row number (int): ' ) ) # Y
                self.superCubeCenterCol = int( input( 'Enter the center column number: ' ) ) # X
                self.voidRatio = float(input('Enter measured void ratio: '))
                self.tiffFileLocation = input( 'Enter the location of tiff files: ' )
                allSampleDetailsObtained = True
            else:
                print( 'The option is incorrect, please re-enter...' ) 
        
    
    # Just for laughs    
    def whatHo(self):
        print("Hello, Sir. How may I be of assistance?")

    def readMeThePoem(self):
        print('Here you go sir...\n')
        print('If you can dream—and not make dreams your master;')   
        print('If you can think—and not make thoughts your aim;')   
        print('If you can meet with Triumph and Disaster')
        print('And treat those two impostors just the same;') 
    