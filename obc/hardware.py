#from machine import Pin
import umatrix

class Magnetometer (object):
	Accuracy = 0
	Resolution = 3.45
	Range = 0
	Readings = []
	I2CBus = None

	def __init__ (self, testRange):
		self.Range = testRange

	def getReadings (self):
		Readings = inputPin.read()
		# TODO: add transformations so all readings line up


class Magnetorquer (object):
	ConductorGuage = 0
	ConductorLength	= 0
	Restivity = 0
	Position = []
	Orientation = []
