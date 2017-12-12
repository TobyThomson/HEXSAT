from collections import namedtuple
from sgp4.propagation import getgravconst

EarthGravity = namedtuple(
    'EarthGravity',
    'tumin mu radiusearthkm xke j2 j3 j4 j3oj2',
    )

wgs84 = EarthGravity(*getgravconst('wgs84'))