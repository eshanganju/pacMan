# -*- coding: utf-8 -*-

'''
	Purpose is to 
		1. Get the REV size that can be used for analysi of all the data
		2. Get particle size parameter that 

'''

import classes.Aggregate as Aggregate
import classes.Reader as Reader
import classes.Filter as Filter
import classes.Segment as Segment
import classes.Measure as Measure
import classes.LemmeC as LemmeC
import classes.Jeeves as Jeeves

# Objects
r = Reader.Reader()
f = Filter.Filter()
s = Segment.Segment()
m = Measure.Measure()
l = LemmeC.LemmeC()
j = Jeeves.Jeeves()

j.checkFolderLocations()
j.checkSampleDetails()
j.readDataFiles()

gliStack = r.readTiffSequence( j.tiffFileLocation, 
							   j.superCubeCenterSlice,
							   j.superCubeCenterRow,
							   j.superCubeCenterCol, 
							   j.superCubeEdgeLengthforREVAnalysis, 
							   j.superCubeCalib )

superCube = Aggregate.Aggregate( j.sandName, 
								 gliStack, 
								 j.superCubeCalib, 
								 16, 
								 j.voidRatio, 
								 True, 
								 j.outputFilesLocation )

f.filterSample( superCube )    

s.binarizeAccordingToDensity( superCube )
s.obtainEuclidDistanceMap( superCube )
s.obtainLocalHMaximaMarkers( superCube , 4)
s.obtainTopoWatershed( superCube )
s.fixErrorsInSegmentation( superCube )
s.resetParticleList( superCube )

m.measureParticleSizeDistribution ( superCube )
m.revSizeAnalysis( superCube )

'''
Correct the void ratrio calculation in binarization otsu
Watershed algorithm - write new code - test in batch. take files for output
Correction of segmentation issues

'''

