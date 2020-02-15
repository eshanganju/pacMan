'''

'''

from classes import PAC as pac

pac = pac.ParticleAnalysisCode()

pac.checkfolderLocations()
pac.checkSampleDetails()
pac.collectDataFiles()

pac.filterData()
pac.segmentData()
pac.measureParticleSizeParams() 
pac.performREVSizeAnalysis()