import tifffile as tf
import time


from pac import Project as p1
from pac import Project_M as p2


image=tf.imread('/home/chawlahpc2adm/Desktop/Daniel/RAMpowder/Output/Retest/Retest-RAM2-Al_643x643x972_h45_r7-noEdgeCorrectedLabelMap-ManualCorrectionmissingLabelsCorrected.tif')

startTimeP1 = time.time()

p1.getTiffstacks( image, 
					saveData=True, 
					fileName='AA7050-RAM2', 
					outputDir = '/home/chawlahpc2adm/Desktop/Daniel/RAMpowder/Output/ClassCheckSTLs/ProjectionsGoSlow/', 
					crop=True, 
					size=250, 
					normalize=True)

endTimeP1 = time.time()

imageStack = p2.getTiffStacks2(image,
								saveData=True, 
								fileName='AA7050-RAM2', 
								outputDir = '/home/chawlahpc2adm/Desktop/Daniel/RAMpowder/Output/ClassCheckSTLs/Scan2Projections/', 
								size=500, 
								normalize=True )

endTimeP2 = time.time()

timeP1 = endTimeP1 - startTimeP1
timeP2 = endTimeP2 - endTimeP1

print( 'Time non parallel: ' + str( timeP1//60) + ' mins' )
print( 'Time parallel: ' + str( timeP2//60 ) + 'mins' )
