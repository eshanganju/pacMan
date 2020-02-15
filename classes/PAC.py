'''

'''

# General imports
import numpy as np
import skimage.external.tifffile as tiffy

# Special imports
from classes import Tangibles as Tangibles
from classes import Filter as Filter
from classes import Segment as Segment
from classes import Measure as Measure
from classes import Reader as Reader

# %% Class
class PAC:

    def __init__(self):
        self.inputFilesLocation = ''
        self.outputFilesLocation = ''
        self.aggregateList = [] 
        print( 'Particle Analysis Code activated!' )
             
    def checkFolderLocations(self):      
        defaultInputFilesLocation = 'C:/Users/eganj/gitHub/pacInput/'
        defaultOutputFilesLocation = 'C:/Users/eganj/gitHub/pacOutput/'       
        
        # Check default input: 
        inputLocationIsCorrect = input( '\n\nBy default, input files are supposed to be at: \n' + defaultInputFilesLocation + '\nIs this correct? ["y"/n]: ' )        
        if inputLocationIsCorrect.lower() == 'n':
            self.inputFilesLocation = input( '\nEnter new location of input files: ' )                
        else:
            self.inputFilesLocation = defaultInputFilesLocation   

        # Check default output: 
        outputLocationIsCorrect = input( '\n\nBy default, output files are supposed to be at: \n' + defaultOutputFilesLocation + '\nIs this correct? ["y"/n]: ' )                 
        if outputLocationIsCorrect.lower() == 'n':
            self.outputFilesLocation = input('\nEnter new location of output files: ')      
        else:
            self.outputFilesLocation = defaultOutputFilesLocation

        print('\nFolder locations established...')
     
    def checkSampleDetails(self):
        print( '\nSample default calibration is 0.01193 mm/px...' )
        changeSampleCalib = input( 'Do you want to change this [y/"n"]?: ' )
        
        if changeSampleCalib.lower() == 'y':
            self.sampleCalib = float( input( 'Enter calibration (mm/px): ' ) )       
        else:
            self.superCubeCalib = 0.01193 
                
        # Supercube size details: 
        print( '\nSample default edge length is 5.5 mm' )
        changeSuperCubeEdgeLength = input( 'Do you want to change this [y/"n"]?: ' )      
        
        if changeSuperCubeEdgeLength.lower() == 'y':
            self.superCubeEdgeLengthforREVAnalysis = float( input( 'Enter superCube edge length (mm): ' ) )
        else:
            self.superCubeEdgeLengthforREVAnalysis = 5.5
                
        # Sand sample details: 
        allSampleDetailsObtained = False       
        while allSampleDetailsObtained == False:
            sand = input( '\nEnter sand name (ogf, 2qr, otc or new): ' )
            if sand.lower() == 'ogf':
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
                self.sandName = input('\nEnter sand name: ')
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
        
    def readData( self ):
        '''
        Creates aggregates
        '''

    def filterData( self ):
        '''
        Filters data using the aggregate
        '''

    def segmentData( self ):
        '''
        Segments data following one of the following: 
            EDT-WS
            Level set
            Random walker      

        '''

    def measureParticleSizeParam( self ):
        '''
        Obtain particle size params of the entire sample
        '''

    def performREVsizeAnalysis( self ):
        '''
        Cycle through different size to obtain REV size
        '''

