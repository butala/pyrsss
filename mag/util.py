import math

"""
The equations are found in the 2010 World Magnetic Model 2015
Report. See (19) on page 11.

http://www.ngdc.noaa.gov/geomag/WMM/data/WMM2015/WMM2015_Report.pdf
"""

def mag_inc(x, y, z):
    """
    Given *x* (north intensity), *y* (east intensity), and *z*
    (vertical intensity) all in [nT], return the magnetic inclincation
    angle [deg].
    """
    h = math.sqrt(x**2 + y**2)
    return math.degrees(math.atan2(z, h))


def mag_dec(x, y, z):
    """
    Given *x*, *y*, *z* as defined in :fun:`mag_inc`, return the
    magnetic declination angle [deg].
    """
    return math.degrees(math.atan2(x, y))
