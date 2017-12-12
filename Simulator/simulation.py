from unittest.mock import Mock
# TODO: use skyfeild (its pure python!!!)
#from pyorbital.orbital import Orbital
import sys
from datetime import datetime, timedelta
from ADC import ADC
import random
import csv
# TODO: use geomag (pure python too!!!!!)
#from pyigrf12 import runigrf12

sys.modules['machine'] = Mock()

from Hardware import Magnetometer, Magnetorquer

class SimulationMagnetometer (Magnetometer):
    Variance = 1
    OrbitPropogorgator = None

    def __init__ (self, variance, orbitPropogorgator, **kargs):
        super().__init__(**kargs)

        self.Variance = variance
        self.OrbitPropogorgator = orbitPropogorgator


    def getReadings(self):
        # TODO: find a neat way of handing this info to the simulation magnetometer (init is probably the way)
        self.orbitPropogorgator.UpdateTime(additionalSeconds=2)

        position = self.orbitPropogorgator.GetPosition()

        latitude = position[1]
        longitude = position[0]
        altitude = position[2]

        #Gets the B-feild value at the simulated position
        igrfValue = runigrf12(simulationTime, 0, 1, altitude, latitude, longitude)
        #Sorts pyigrf12 weirdness
        igrfValues = [x for y in igrfValues for x in y][:3]
        #Adds gaussain white noise to each reading
        igrfValues = [random.gauss(x, variance) for x in igrfValues]
        #nT -> mG
        igrfValues = [x* 10**-2 for x in igrfValues]
        #Rounds each reading according to the magnometer's resolution
        self.Readings = [round(x / magnetometer.Resolution) * magnetometer.Resolution for x in igrfValues]

        return self.Readings

class SimulationMagnetorquer (Magnetorquer):
    def Pulse(self, current, duration):
        # TODO: log to csv file
        print('Just pulsed!')

class OrbitPropogorgator:
    def __init__ (self):
        self.SimulationTime = datetime.utcnow()
        self.Orbit = Orbital("NOAA 18")

    def UpdateTime (self, additionalSeconds=0, additionalMinutes=0):
        self.SimulationTime += timedelta(seconds=additionalSeconds, minutes=additionalMinutes)

    def GetPosition (self):
        position = orbit.get_lonlatalt(self.SimulationTime)
        return position

#***************************SETUP*********************************
mOrbitPropgator = OrbitPropogorgator()

SimulationDurationHours = 6
MagnetometerReadingVariance = 500

SimulationMagnetometers = []
SimulationMagnetorquers = []

for x in range(3):
    SimulationMagnetometers.append(SimulationMagnetometer(testRange=500, variance=200, orbitPropogorgator=mOrbitPropgator))
    SimulationMagnetorquers.append(SimulationMagnetorquer())

SimulationADC = ADC(SimulationMagnetometers, SimulationMagnetorquers)

for minute in range(SimulationDurationHours * 60):
    mOrbitPropgator.UpdateTime(additionalMinutes=1)
    SimulationADC.Update()
