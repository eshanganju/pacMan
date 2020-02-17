'''
PAC main
'''

from classes import PAC as PAC

pac = PAC.PAC()

pac.readNewData()
pac.filterData()
pac.segmentData()