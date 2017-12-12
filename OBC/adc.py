import umatrix
import ulinalg

class ADC:
    isTumbling = True

    Magnetometers = []
    Magnetorquers = []

    BDotGain = 10**5
    MaxCurrent = 0.2
    SampleTime = 2.0
    N = 3.5
    A = 0.001

    def __init__ (self, magnetometers, magnetorquers):
        self.Magnetometers = magnetometers
        self.Magnetorquers = magnetorquers

    def GetMagnetometerReadings (self):
        for x in range(len(self.Magnetometers)):
            readings = self.Magnetometers[x].getReadings()
            adverageReadings[x,:] = umatrix.matrix([readings])

        #adveraging stuff

        return adverageReadings

    def BDotControl (self):
        #issue is that the second readings are the same
        firstReadings = self.GetMagnetometerReadings()
        print(firstReadings)
        #sleep(self.SampleTime)
        secondReadings = self.GetMagnetometerReadings()

        bDot = (firstReadings - secondReadings).reciprocal() * self.SampleTime

        print(bDot)
        m = np.multiply(-self.BDotGain, bDot)
        print(m)
        i = np.divide(m, (self.N, self.A))

        i = np.clip(i, -self.MaxCurrent, self.MaxCurrent)

        for x in range(len(self.Magnetorquers)):
            self.Magnetorquers[x].Pulse(i[x], self.PulseLength)

    def Update (self):
        if (self.isTumbling):
            self.BDotControl()
