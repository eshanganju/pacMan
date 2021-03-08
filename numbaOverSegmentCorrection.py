"""This code is to implement numba in the oversegmentation correction
The objective is to write one self-contained code that uses only numpy and numba
"""

import numba as nb
import numpy as np

def fixOverSegmentation( labMap,areaLimit=0,radiusRatioLimit=0,saveLog=True,
                         saveImg=True,sampleName='',outputDir='' ):
    """This function corrects oversegmentation by merging labels based on either a
    contact area limit or a radius ratio limit

    Parameters
    ----------
    labMap : 
    areaLimit :
    radiusRatioLimit :
    saveLog : 
    saveImg : True
    sampleName : 
    outputDir :  
    
    Return
    ------
    correctedLabelMap : 

    """
