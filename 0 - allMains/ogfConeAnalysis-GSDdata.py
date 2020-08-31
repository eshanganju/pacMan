import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt
from classes import Measure

inputDataLocation = '/home/eg/codes/pacInput/zaManualSegmentation2020-03-04/editedLabelFiles/'

searchString = inputDataLocation + '/*csv'

readAllinFolder = True
if readAllinFolder == True:
    fileList = glob.glob(searchString)
else:
    fileList = ['/home/eg/codes/pacInput/zaManualSegmentation2020-03-04/editedLabelFiles/ogf-25kPa-tip-spamLabelledMap-edited.csv']

for fileNum in range(0,len(fileList)):
    indf = pd.read_csv(fileList[fileNum])
    indf.rename(columns={indf.columns[0] : 'ptclNum', \
                         indf.columns[1] : 'vol', \
                         indf.columns[2] : 'eqsp', \
                         indf.columns[3] : 'caMax', \
                         indf.columns[4] : 'caMed', \
                         indf.columns[5] : 'caMin', \
                         indf.columns[6] : 'na1', \
                         indf.columns[7] : 'na2'}, inplace=True)
    
    outdf = Measure.computeGrainSizeDistribution(indf)
    fileName = (fileList[fileNum].split('/')[-1]).split('.')[0]
    completeFileName = inputDataLocation + fileName + 'GSD.csv'
    outdf.to_csv(completeFileName,index=False)
