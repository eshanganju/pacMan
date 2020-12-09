'''
Code for the analysis of 3D tomo data for samples collected around the cone penetrometer

written by @eg

Steps:
    -Initialize the samples variables
    -Read the files and extract subregions from the 3D tomo data - 6D50 edge length

    -Compute particle sizes of the sand samples - user threshold after initial guess from OTSU
    -Get particel aspect ration from PCA length ratios
    -Get particle size distribution
    -Calculate relative breakage parameter
    -Get contact normals
    -Get fabric tensors
    -Get fabric projection plots

'''

from pac import Reader
from pac import Segment
from pac import Measure
from pac import Plot

scanData = ['2QR_25_top','2QR_25_mid','2QR_25_tip',
            '2QR_50_top','2QR_50_mid','2QR_50_tip',
            '2QR_90_top','2QR_90_mid','2QR_90_tip',
            'OGF_25_top','OGF_25_mid','OGF_25_tip',
            'OGF_50_top','OGF_50_mid','OGF_50_tip',
            'OGF_90_top','OGF_90_mid','OGF_90_tip',
            'OTC_25_top','OTC_25_mid','OTC_25_tip',
            'OTC_50_top','OTC_50_mid','OTC_50_tip',
            'OTC_90_top','OTC_90_mid','OTC_90_tip']


for each scan in scanData:
    scanInputLoc = '/home/eg/codes/pacInput/' + scan + '/'
    subregionZYXListLoc = '/home/eg/codes/pacInput/' + scan + '/' + 'subregionZYXLocs.csv'
    '''
    D50 and Calibfactor code selection code
    '''
    subregionYXArray = np.genfromtxt(subregionZYXListLoc, delimiter=',', skip_header=2)
    subregionZ = int(np.genfromtxt(subregionZYXListLoc, delimiter=',', skip_footer=10))
    numberSubregions = subregionXYList.shape[0]

    analysisOutputLoc = '/home/eg/codes/pacOutput/cone/' + scan + '/'
    imageOutputLoc = '/media/eg/EshanGanju/coneImages/' + scan + '/'

    for currentSubregion in range(0,numberSubregions):
        subregionGLIMap = Reader.tifffileSequence2(folderLocation=scanInputLoc,
                                                   centerZ=subregionZ,
                                                   topLeftY=subregionYXArray[currentSubregion,0],
                                                   topLeftX=subregionYXArray[currentSubregion,1],
                                                   lngt=6*D50,
                                                   calib=calibrationFactor,
                                                   invImg=False)
