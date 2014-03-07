from math import radians
from numpy import pi, arctan2, sin, cos, arcsin
import sys

def gal2eq(l, b):
    # convert first from degrees to radians
    l = radians(l)
    b = radians(b)
    
    # Here we convert from Galactic to Equatorial (all conversions done in radians)
    # North galactic pole (J2000) -- according to Wikipedia
    pole_ra = radians(192.859508)
    pole_dec = radians(27.128336)
    posangle = radians(122.932-90.0)
    ra = arctan2( (cos(b)*cos(l-posangle)), (sin(b)*cos(pole_dec) - cos(b)*sin(pole_dec)*sin(l-posangle)) ) + pole_ra
    dec = arcsin( cos(b)*cos(pole_dec)*sin(l-posangle) + sin(b)*sin(pole_dec) )

    # Convert from radians to degrees for output in degrees
    ra = ra * 180./pi # degrees
    dec = dec * 180./pi # degrees

    return (ra, dec)


