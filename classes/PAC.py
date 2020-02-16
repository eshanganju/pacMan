'''
PAC  class
This is like the controller. All the subsequest classes are accessed through this file. 
This class can do the following: 
    1. 
'''

# General imports
import numpy as np
import skimage.external.tifffile as tiffy

# Special imports
from classes import Aggregate as Aggregate
from classes import Filter as Filter
from classes import Segment as Segment
from classes import Measure as Measure
from classes import Reader as Reader

# %% Class
class PAC:

    def __init__(self):
        self.inputFilesLocation = ''
        self.outputFilesLocation = ''
        self.currentAggregateListIndex = 0
        self.totalNumberOfAggregates = 0
        self.aggregateList = [] 
        print( '\n-------------------------------*' )
        print( 'Particle Analysis Code activated!' )
        print( '-------------------------------*' )

    def readNewData( self ):

        self.checkFolderLocations()
        smplName, btdpth, smplClb, cubEdgLngt, smplCntrZ, smplCntrY, smplCntrX, d50, voidRatio, tifFilFldrLoc = self.getSampleDetails()
        
        # Read GLI data
        r = Reader.Reader()
        dataFormTiffStack = input( '\nIs the data in the form of a tiff stack (y/[n])?: ' )
        if dataFormTiffStack == 'y':
            tifStkFilNam = input( 'Enter name of the tiff stack file: ' )
            tifStkFilAndFldrNam = self.inputFilesLocation + tifStkFilNam
            gliDat = r.readTiffStack( tifStkFilAndFldrNam, smplCntrZ, smplCntrY, smplCntrX, cubEdgLngt, smplClb ) 
        else:
            gliDat = r.readTiffFileSequence( tifFilFldrLoc, smplCntrZ, smplCntrY, smplCntrX, cubEdgLngt, smplClb ) 

        # Create Aggregate
        self.aggregateList[ self.currentAggregateListIndex ] = Aggregate.Aggregate( smplName, btdpth, smplClb, cubEdgLngt, smplCntrZ, smplCntrY, smplCntrX, d50, voidRatio, tifFilFldrLoc, gliDat)
        
        # Update directory
        self.currentAggregateListIndex = self.currentAggregateListIndex + 1
        self.totalNumberOfAggregates = self.totalNumberOfAggregates + 1
        
        # Cleanup
        del r
        del gliDat



    def checkFolderLocations(self):      
        print('\n\nChecking folder locations: ')
        print('---------------------------*')

        defaultInputFilesLocation = 'C:/Users/eganj/gitHub/pacInput/'
        defaultOutputFilesLocation = 'C:/Users/eganj/gitHub/pacOutput/'       
        
        # Check default input: 
        inputLocationIsCorrect = input( '\nBy default, input files are supposed to be at: \n' + defaultInputFilesLocation + ', Is this correct? ([y]/n): ' )        
        if inputLocationIsCorrect.lower() == 'n':
            self.inputFilesLocation = input( '\nEnter new location of input files: ' )                
        else:
            self.inputFilesLocation = defaultInputFilesLocation   

        # Check default output: 
        outputLocationIsCorrect = input( '\n\nBy default, output files are supposed to be at: \n' + defaultOutputFilesLocation + '\nIs this correct? ([y]/n): ' )                 
        if outputLocationIsCorrect.lower() == 'n':
            self.outputFilesLocation = input('\nEnter new location of output files: ')      
        else:
            self.outputFilesLocation = defaultOutputFilesLocation

        print('\nFolder locations established...')
     

    def getSampleDetails(self):
        
        print('\n\nChecking sample details: ')
        print('---------------------------*')

        # Sample name: 
        sampleName = input( '\nEnter the name of the sample: ' )

        # Bitdepth: 
        bitDepth = input( '\nEnter bitdepth of the images (generally 8, 16 or 32): ' )

        # Calibration: 
        print( '\nSample default calibration is 0.01193 mm/px...' )
        changeSampleCalib = input( 'Do you want to change this (y/[n]): ' )
        if changeSampleCalib.lower() == 'y':
            sampleCalib = float( input( 'Enter calibration (mm/px): ' ) )       
        else:
            sampleCalib = 0.01193 
                
        # Size: 
        print( '\nSample default edge length of cube analyzed is 5.5 mm' )
        changeCubeEdgeLength = input( 'Do you want to change this (y/[n]): ' )      
        if changeCubeEdgeLength.lower() == 'y':
            cubeEdgeLength = float( input( 'Enter new edge length (mm): ' ) )
        else:
            cubeEdgeLength = 5.5
        
        # Center location:         
        print('\nChecking center of the slices: ')
        sampleCenterZ = int(input( 'Location of center slice of the sample (px): ' ))
        sampleCenterY = int(input( 'Location of center row of the sample (px): ' ))
        sampleCenterX = int(input( 'Location of center column of the sample (px): ' ))
        
        # Average particle size: 
        d50Known = input( '\nIs the original D50 of sample know ([y]/n): ' )
        if d50Known == 'n':
            d50 = np.nan
        else: 
            d50 = float(input( 'Enter the original D50 of sample (mm): ' ))

        # Void ratio
        voidRatioKnown = input( '\nIs the sample void ratio know ([y]/n): ' )
        if voidRatioKnown == 'n':
            voidRatio = np.nan
        else: 
            voidRatio = float( input( 'Enter the void ratio of sample (decimals): ' ) )

        # Data folder name in pacInput fodler: 
        tiffFileFolderName = input('\nName of data folder in ' + self.inputFilesLocation + ' where files are located: ')
        tiffFileFolderLocation = self.inputFilesLocation + tiffFileFolderName

        return sampleName, bitDepth, sampleCalib, cubeEdgeLength, sampleCenterZ, sampleCenterY, sampleCenterX, d50, voidRatio, tiffFileFolderLocation


    def filterData( self ):
        '''
        Filters data using the aggregate
        '''

    def segmentData( self ):
        '''
        Segments data following one of the following
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

