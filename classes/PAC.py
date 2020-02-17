'''
PAC  class
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
        self.totalNumberOfAggregates = 0
        self.aggregateList = [] 
        print( '\n-------------------------------*' )
        print( 'Particle Analysis Code activated!' )
        print( '-------------------------------*' )

    def readNewData( self ):

        self.checkFolderLocations()
        smplName, btdpth, smplClb, cubEdgLngt, smplCntrZ, smplCntrY, smplCntrX, d50, voidRatio, datFrmTifStk, tifFilFldrLoc = self.getSampleDetails()
        
        # Read GLI data
        r = Reader.Reader()
        
        if datFrmTifStk == True:           
            gliDat = r.readTiffStack( tifFilFldrLoc, smplCntrZ, smplCntrY, smplCntrX, cubEdgLngt, smplClb )
            self.aggregateList.append(Aggregate.Aggregate( smplName, btdpth, smplClb, cubEdgLngt, smplCntrZ, smplCntrY, smplCntrX, d50, voidRatio, datFrmTifStk, tifFilFldrLoc, gliDat ))
        else:
            gliDat = r.readTiffFileSequence( tifFilFldrLoc, smplCntrZ, smplCntrY, smplCntrX, cubEdgLngt, smplClb ) 
            self.aggregateList.append(Aggregate.Aggregate( smplName, btdpth, smplClb, cubEdgLngt, smplCntrZ, smplCntrY, smplCntrX, d50, voidRatio, datFrmTifStk, tifFilFldrLoc, gliDat ))

        # Update directory
        self.totalNumberOfAggregates = self.totalNumberOfAggregates + 1
        
        # Create Aggregate
        createGLIFile = input('Do you want to save a copy of the GLI as a stack in the output foder (y/[n])?')
        if createGLIFile.lower() == 'y':
            unfilteredImageName = str(self.outputFilesLocation) + str( smplName ) + '.tiff'
            tiffy.imsave( unfilteredImageName, self.aggregateList[ self.totalNumberOfAggregates - 1 ].greyLevelMap )
            print( "\nUnfiltered GLI file saved as " + str( smplName ) + '.tiff' )         
        else:
            print('Ok, File not saved')
        
        print('Total number of aggregates is %d\n' % self.totalNumberOfAggregates)
        print('Data reading complete------------*\n\n')
        
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
        outputLocationIsCorrect = input( '\n\nBy default, output files will be stored in: \n' + defaultOutputFilesLocation + '\nIs this correct? ([y]/n): ' )                 
        if outputLocationIsCorrect.lower() == 'n':
            self.outputFilesLocation = input('\nEnter new location of output files: ')      
        else:
            self.outputFilesLocation = defaultOutputFilesLocation

        print('\nFolder locations established...')

    def getSampleDetails(self):
        
        print('\n\nChecking sample details: ')
        print('---------------------------*')

        # Sample name: 
        useDefaultSampleName = input('\nIs this OTC-0N sample? ([y]/n):')
        if useDefaultSampleName == 'n':
            sampleName = input( 'Enter the name of the sample: ' )
        else:
            sampleName = 'otc'

        # Bitdepth: 
        useDefaultBitDepth = input('\nUse default bit depth of 16? ([y]/n)')
        if useDefaultBitDepth == 'n':
            bitDepth = int( input( 'Enter bitdepth of the images (generally 8, 16 or 32): ' ) )
        else:
            bitDepth = 16

        # Calibration: 
        useDefaultSampleCalib = input( '\nUse default calibration of 0.01193 mm/px? ([y]/n): ' )
        if useDefaultSampleCalib.lower() == 'n':
            sampleCalib = float( input( 'Enter calibration (mm/px): ' ) )       
        else:
            sampleCalib = 0.01193 
                
        # Size: 
        useDefaultCubeEdgeLength = input( '\nUse default cube edge lenght of 5.5 mm? ([y]/n): ' )      
        if useDefaultCubeEdgeLength.lower() == 'n':
            cubeEdgeLength = float( input( 'Enter new edge length (mm): ' ) )
        else:
            cubeEdgeLength = 5.5
        
        # Center location:         
        useDefaultOtcCenterPoints = input('\nUse default center of the slices? [Z = 513][Y = 430][X = 490] ([y]/n)')
        if useDefaultOtcCenterPoints == 'n':
            sampleCenterZ = int(input( 'Location of center slice of the sample (px): ' ))
            sampleCenterY = int(input( 'Location of center row of the sample (px): ' ))
            sampleCenterX = int(input( 'Location of center column of the sample (px): ' ))
        else: 
            sampleCenterZ = 513
            sampleCenterY = 430
            sampleCenterX = 490
        
        # Average particle size: 
        d50Known = input( '\nIs the original D50 of sample know ([y]/n): ' )
        if d50Known == 'n':
            d50 = np.nan
        else: 
            useDefaultOTCd50 = input('Use default OTC D50 [0.72 mm] ([y]/n):')
            if useDefaultOTCd50 == 'n':
                d50 = float(input( 'Enter the original D50 of sample (mm): ' ))
            else: d50 = 0.72

        # Void ratio
        voidRatioKnown = input( '\nIs the current sample void ratio know ([y]/n): ' )
        if voidRatioKnown == 'n':
            voidRatio = np.nan
        else: 
            useDefaultOTCvoidRatio = input('Use default OTC void ratio [0.5404] ([y]/n):')     
            if useDefaultOTCvoidRatio == 'n':       
                voidRatio = float( input( 'Enter the current void ratio of the sample (decimals): ' ) )
            else: voidRatio = 0.5404

        # Data folder name in pacInput folder: 
        dataFromTiffStack = input( '\nIs the data in the form of a tiff stakc (y/[n]): ' )
        if dataFromTiffStack.lower() == 'y':
            dataFromTiffStack = True
            print('Make sure the tiff stack is stored in the pacInput folder.')
            tifStkFilNam = input( 'Enter name of the tiff stack file: ' )
            tiffFileFolderLocation = self.inputFilesLocation + tifStkFilNam + '.tiff'
        else: 
            dataFromTiffStack = False
            useDefault0NFolderForOTCSample = input('Use the default folder [OTC-0N]? ([y]/n): ')
            
            if useDefault0NFolderForOTCSample == 'n':
                tiffFileFolderName = input('\nName of data folder in ' + self.inputFilesLocation + ' where tiff file sequences are located: ')
            
            else: 
                tiffFileFolderName = 'OTC-0N'

            tiffFileFolderLocation = self.inputFilesLocation + tiffFileFolderName

        return sampleName, bitDepth, sampleCalib, cubeEdgeLength, sampleCenterZ, sampleCenterY, sampleCenterX, d50, voidRatio, dataFromTiffStack, tiffFileFolderLocation

    def filterData( self ):
        '''
        Filters data using the aggregate
        '''
        print( '\n\nStarting filtration module' )
        print( '----------------------------*\n' )

        filteredDataAvailable = input( 'Is filtered data available? (y/[n]): ' )
        
        if filteredDataAvailable == 'y':
            print('\nDefault filter data location is: ' + self.inputFilesLocation)
            defaultFilteredDataLocationCorrect = input( 'Is this correct? ([y]/n): ' )
            
            if defaultFilteredDataLocationCorrect == 'n':
                newFilteredDataLocation = input( 'Enter filtered tiff stack location: ' )
                newFilteredFileName = input('Enter the name of the file located in ' + newFilteredDataLocation + ' :')
                print('Reading GLI data from new location...')
                newFilteredFileLocation = newFilteredDataLocation + newFilteredFileName
                self.aggregateList[0].filteredGreyLevelMap = tiffy.imread( newFilteredFileLocation )
            
            else:
                filteredFileName = input( 'Enter the name of the file located in ' + self.inputFilesLocation )
                filteredFileLocation = self.inputFilesLocation + filteredFileName
                print('Reading GLI data from default location')
                self.aggregateList[0].filteredGreyLevelMap = tiffy.imread( filteredFileLocation )
            
            self.aggregateList[0].imageNoise = abs( self.aggregateList[0].greyLevelMap - self.aggregateList[0].filteredGreyLevelMap )
            
            # Keeping records
            finalFilteredFileLocation = self.outputFilesLocation + self.aggregateList[0].sampleName + '-filtered.tiff'       
            tiffy.imsave(finalFilteredFileLocation, self.aggregateList[0].filteredGreyLevelMap)

            finalNoiseFileLocation =  self.outputFilesLocation + self.aggregateList[0].sampleName +  '-noise.tiff'            
            tiffy.imsave(finalNoiseFileLocation, self.aggregateList[0].imageNoise)    
        
        else:
            f = Filter.Filter()
            fileName = self.aggregateList[0].sampleName
            gliMap = self.aggregateList[0].greyLevelMap
            gliMaxVal = self.aggregateList[0].gliMax
            outputDirectory = self.outputFilesLocation
            
            self.aggregateList[0].filteredGreyLevelMap = f.filterDenoiseNlm(fileName, gliMap, gliMaxVal, outputDirectory)
            self.aggregateList[0].imageNoise = abs( self.aggregateList[0].greyLevelMap - self.aggregateList[0].filteredGreyLevelMap )
            
            # Keeping records
            filteredFileName = outputDirectory + fileName + '-filtered.tiff'
            noiseFileName =  outputDirectory + fileName +  '-noise.tiff'
            tiffy.imsave(filteredFileName, self.aggregateList[0].filteredGreyLevelMap)
            tiffy.imsave(noiseFileName, self.aggregateList[0].imageNoise)

    def segmentData( self ):
        '''
        Segments data following one of the following
            EDT-WS
            Level set
            Random walker      
        '''
        print( '\n\nStarting segmentation module' )
        print( '----------------------------*\n' )       

        s = Segment.Segment()        
        
        labelledDataAvailable = input( 'Is labelled data available? (y/[n]):' )
        if labelledDataAvailable = 'y':
            print( 'Default labelled data location is: ' + self.inputFilesLocation )
            defaultLabelledDataLocationCorrect = input( 'Is this correct? ([y]/n): ' )

            if defaultLabelledDataLocationCorrect == 'n':
                newLabelledDataLocation = input('Enter labelled tiff stack location: ')
                newLabelledFileName = input('Enter name of the file located in ' + newLabelledDataLocation + ':')
                print('Reading labelled data from new location...')
                newLabelledFileLocation = newLabelledDataLocation + newLabelledFileName
                self.aggregateList[0].labelledMap = tiffy.imread( newLabelledFileLocation )
            
            else:
                labelledFileName = input( 'Enter the name of the file located in ' + self.inputFilesLocation )
                labelledFileLocation = self.inputFilesLocation + labelledFileName
                print('Reading labelled data from location...')
                self.aggregateList[0].labelledMap = tiffy.imread( labelledFileLocation )

            # Keeping record
            finalLabelledFileLocation = self.outputFilesLocation + self.aggregateList[0].sampleName + '-labelled.tiff' 
            tiffy.imsave(finalLabelledFileLocation, self.aggregateList[0].labelledMap)
        else:
            s=Segment.Segment()

            print('\nSegmentation method: ')
            segmentationMethod = int( input( 'Choose: ([1]) EDT-WS, (2) LS-WS, (3) RW (choose 1, 2 or 3): ') )

            # Random Walker watershed
            if segmentationMethod == 3:
                print('\nStarting Random Walker Segmentation...')

            # Level set segmentation
            elif segmentationMethod == 2:
                print('\nStarting Level Set Segmentation...')

            # Euclidian distance map watershed segmentation
            else:
                print('\nStarting EDT - Watershed segmentation...')
                binMap, voidRatioCT, edMap, mrkrMap, lblCorrectionMethod, correctedLabelledMap = s.performEDTWS( self.aggregateList[0].filteredGreyLevelMap, 
                                                                                                                 self.aggregateList[0].currentvoidRatioMeasured, 
                                                                                                                 self.outputFilesLocation,
                                                                                                                 self.sampleName ) 
                self.aggregateList[0].binaryMap = binMap
                self.aggregateList[0].currentvoidRatioFromCT = voidRatioCT
                self.aggregateList[0].euclidDistanceMap = edMap
                self.aggregateList[0].edPeakMarkers = mrkrMap
                self.aggregateList[0].labelledMap = correctedLabelledMap

                del binMap, voidRatioCT, edMap, mrkrMap, correctedLabelledMap

                finalBinaryFileLocation = self.outputFilesLocation + self.aggregateList[0].sampleName + '-binary.tiff'
                finalEDMLocation = self.outputFilesLocation + self.aggregateList[0].sampleName + '-edm.tiff'
                finalEDPeakMarkersLocation = self.outputFilesLocation + self.aggregateList[0].sampleName + '-edmPeaks.tiff'
                finalLabelledMapLocation = self.outputFilesLocation + self.aggregateList[0].sampleName + '-labelledMap.tiff'

                tiffy.imsave(finalBinaryFileLocation, self.aggregateList[0].binaryMap)
                tiffy.imsave(finalEDMLocation, self.aggregateList[0].euclidDistanceMap)
                tiffy.imsave(finalEDPeakMarkersLocation, self.aggregateList[0].edPeakMarkers)
                tiffy.imsave(finalLabelledMapLocation, self.aggregateList[0].labelledMap)

    def measureParticleSizeParam( self ):
        '''
        Obtain particle size params of the entire sample
        '''

    def performREVsizeAnalysis( self ):
        '''
        Cycle through different size to obtain REV size
        '''

    def plotCTSlices(self, imageDataAs3DNumpyArray):
        '''
        Plotting slices
        '''



