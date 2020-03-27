'''
Particle class
'''
import numpy as np

class Particle:

    def __init__(self, num, numVox, locData):
        self.index = num
        self.volume = numVox
        self.locationData = locData
        self.covarianceMatrix = np.zeros((3,3))
        self.eigenValue = np.zeros(3)
        self.eigenVector = np.zeros((3,3))
        self.meanZ = 0
        self.meanY = 0
        self.meanX = 0
        self.centeredLocationData = np.zeros_like(self.locationData)
        self.rotatedCenteredLocationData = np.zeros_like(self.locationData)
        self.feretMax = 0
        self.feretMin = 0
        self.feretMed = 0
        self.equivalentSphereDiameter = 0
        print("Particle No.",self.index,"created")
