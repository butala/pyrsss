import logging
from collections import defaultdict

import pandas as PD

from ..util.fortran import descriptor2colspecs


"""
From http://omniweb.gsfc.nasa.gov/html/HROdocum.html, retrieved
2017-03-11.

YYYY DDD - YYYY DDD  MM/DD    Parameters
---------------------------
1981 001 - 1994 365 (12/31)   IMF (from IMP8 only)
1981 001 - 1994 365 (12/31)   Plasma (from IMP8 only)
---------------------------------------

1995 001 - 2017 049 (02/18)   IMF and Plasma (from IMP8, ACE, Wind)

1981 001 - 1988 182 (06/30)   Final AE, AL, AU indexes
1988 183 - 1989 365 (12/31)   No AE, AL,AU-indexes (except March of 1989)
1990 001 - 2016 366 (12/31)   Provisional AE, AL, AU indexes

1981 001 - 2017 059 (02/18)   Provisional SYM/D, SYS/H, ASYM/D, ASYS/H indexes

1981 001 - 2014 365 (12/31)   PCN index  (Preliminary)

1986 001 - 2017 059 (02/18)   Fluxes from Goes, for 5-min res. only
---------------------------------------------------------------------

Time span of the Spacecraft Specific Data shifted to bow shock nose
Geotail:  1995-03-15 - 2006-12-31(365)
IMP-8:    1973-11-04 - 2000-06-09(161)
ACE:      1998-02-05 - 2016-12-25(360)
Wind:     1995-01-01 - 2017-02-18(049)


4b. High-Resolution OMNI data set

The common format for the 1-min and 5-min OMNI data sets is

Year			        I4	      1995 ... 2006
Day			        I4	1 ... 365 or 366
Hour			        I3	0 ... 23
Minute			        I3	0 ... 59 at start of average
ID for IMF spacecraft	        I3	See  footnote D below
ID for SW Plasma spacecraft	I3	See  footnote D below
# of points in IMF averages	I4
# of points in Plasma averages	I4
Percent interp		        I4	See  footnote A above
Timeshift, sec		        I7
RMS, Timeshift		        I7
RMS, Phase front normal	        F6.2	See Footnotes E, F below
Time btwn observations, sec	I7	DBOT1, See  footnote C above
Field magnitude average, nT	F8.2
Bx, nT (GSE, GSM)		F8.2
By, nT (GSE)		        F8.2
Bz, nT (GSE)		        F8.2
By, nT (GSM)	                F8.2	Determined from post-shift GSE components
Bz, nT (GSM)	                F8.2	Determined from post-shift GSE components
RMS SD B scalar, nT	        F8.2
RMS SD field vector, nT	        F8.2	See  footnote E below
Flow speed, km/s		F8.1
Vx Velocity, km/s, GSE	        F8.1
Vy Velocity, km/s, GSE	        F8.1
Vz Velocity, km/s, GSE	        F8.1
Proton Density, n/cc		F7.2
Temperature, K		        F9.0
Flow pressure, nPa		F6.2	See  footnote G below
Electric field, mV/m		F7.2	See  footnote G below
Plasma beta		        F7.2	See  footnote G below
Alfven mach number		F6.1	See  footnote G below
X(s/c), GSE, Re		        F8.2
Y(s/c), GSE, Re		        F8.2
Z(s/c), GSE, Re		        F8.2
BSN location, Xgse, Re	        F8.2	BSN = bow shock nose
BSN location, Ygse, Re	        F8.2
BSN location, Zgse, Re 	        F8.2

AE-index, nT                    I6      See World Data Center for Geomagnetism, Kyoto
AL-index, nT                    I6      See World Data Center for Geomagnetism, Kyoto
AU-index, nT                    I6      See World Data Center for Geomagnetism, Kyoto
SYM/D index, nT                 I6      See World Data Center for Geomagnetism, Kyoto
SYM/H index, nT                 I6      See World Data Center for Geomagnetism, Kyoto
ASY/D index, nT                 I6      See World Data Center for Geomagnetism, Kyoto
ASY/H index, nT                 I6      See World Data Center for Geomagnetism, Kyoto
PC(N) index,                    F7.2    See World Data Center for Geomagnetism, Copenhagen
Magnetosonic mach number        F5.1    See  footnote G below
-------
Proton flux (>10 MeV)           F9.2    In 5-min OMNI, but not in 1-min OMNI
Proton flux (>30 MeV)           F9.2    In 5-min OMNI, but not in 1-min OMNI
Proton flux (>60 MeV)           F9.2    In 5-min OMNI, but not in 1-min OMNI


The data may be read with the format statement
1-min: (2I4,4I3,3I4,2I7,F6.2,I7, 8F8.2,4F8.1,F7.2,F9.0,F6.2,2F7.2,F6.1,6F8.2,7I6,F7.2,F5.1)
5-min: (2I4,4I3,3I4,2I7,F6.2,I7, 8F8.2,4F8.1,F7.2,F9.0,F6.2,2F7.2,F6.1,6F8.2,7I6,F7.2,F5.1,3F9.2)

Note that for missing data, fill values consisting of a blank followed
by 9's which together constitute the Ix or Fx.y format are used.

Footnote D:

The following spacecraft ID's are used:
	ACE	71
	Geotail	60
	IMP 8	50
	Wind	51

Footnote E:

Note that standard deviations for the minute-averaged phase front
normal and magnetic field vectors are given as the square roots of the
sum of squares of the standard deviations in the component averages.
For the magnetic field vectors only, the component averages are given
in the records but not their individual standard deviations.  1-min
averaged phase front normal directions are given in the
spacecraft-specific data sets but not in the high resolution OMNI data
set.

Footnote F:

There are no phase front normal standard deviations in the 5-min
records.  This word has fill (99.99) for such records.

Footnote G:

Derived parameters are obtained from the following equations.

Flow pressure = (2*10**-6)*Np*Vp**2 nPa (Np in cm**-3, Vp in km/s,
subscript "p" for "proton")

Electric field = -V(km/s) * Bz (nT; GSM) * 10**-3

Plasma beta = [(T*4.16/10**5) + 5.34] * Np / B**2 (B in nT) (Note that
very low |B| values (<~ 0.3 nT) encountered rarely in high resolution
data can drive plasma beta values to above 1000.  In high resolution
OMNI, there were about 20 such minutes encountered in ~12 years.  We
have assigned the value 998.0 to plasma beta in such cases.  Correct
values of T, Np and B are available in the records for recomputation
of plasma beta values.)

Alfven Mach number = (V * Np**0.5) / (20 * B)

Magnetosonic Mach Number = V/Magnetosonic_speed Magnetosonic speed =
[(sound speed)**2 + (Alfv speed)**2]**0.5 The Alfven speed = 20. * B /
N**0.5 The sound speed = 0.12 * [T + 1.28*10**5]**0.5

For details on these, see https://omniweb.gsfc.nasa.gov/ftpbrowser/bow_derivation.html
"""


"""Fortran string description of the fixed format file."""
COLSPECS_1MIN_STR = '2I4,4I3,3I4,2I7,F6.2,I7, 8F8.2,4F8.1,F7.2,F9.0,F6.2,2F7.2,F6.1,6F8.2,7I6,F7.2,F5.1'


"""
Column specifications suitable for pandas fixed format reader.
"""
COLSPECS = descriptor2colspecs(COLSPECS_1MIN_STR)


"""
Column names.
"""
NAMES = ['year',                 # 1
         'day',
         'hour',
         'minute',
         'IMF_id',
         'SW_id',
         'N_IMF',
         'N_plasma',
         'percent_interp',
         'timeshift',            # 10
         'RMS_timeshift',
         'RMS_phase_front_normal',
         'delta_obs',
         'B_avg',
         'Bx',
         'By_GSE',
         'Bz_GSE',
         'By_GSM',
         'Bz_GSM',
         'RMS_SD_B_scalar',      # 20
         'RMS_SD_B_vector',
         'flow_speed',
         'Vx',
         'Vy',
         'Vz',
         'proton_density',
         'temperature',
         'flow_pressure',
         'electric_field',
         'plasma_beta',          # 30
         'alfven_mach_number',
         'X_GSE',
         'Y_GSE',
         'Z_GSE',
         'BSN_X_GSE',
         'BSN_Y_GSE',
         'BSN_Z_GSE',
         'AE',
         'AL',
         'AU',                   # 40
         'SYM_D',
         'SYM_H',
         'ASY_D',
         'ASY_H',
         'PC_N',
         'magnetosonic_mach_number']  # 46


"""
Example line containing acceptable NaN values. Columns beginning
with X can never have a NaN.
"""
NAN_LINE = 'XXXX  XX XX XX 99 99 999 999 999 999999 999999 99.99 999999 9999.99 9999.99 9999.99 9999.99 9999.99 9999.99 9999.99 9999.99 99999.9 99999.9 99999.9 99999.9 999.99 9999999. 99.99 999.99 999.99 999.9 9999.99 9999.99 9999.99 9999.99 9999.99 9999.99 99999 99999 99999     X     X     X     X 999.99 99.9'


"""
Mapping between column names and acceptable NaN values.
"""
NA_VALUES = defaultdict(list)
for name, tok in zip(NAMES, NAN_LINE.split()):
    if tok.startswith('9'):
        NA_VALUES[name] = [tok]


def parse(omni_fname,
          colspecs=COLSPECS,
          names=NAMES,
          na_values=NA_VALUES):
    """
    Parse the OMNI data record *omni_fname* and return a
    :class:`DataFrame`. To parse, use the fixed columns *colspecs*,
    the column identifiers *names*, and acceptable NaN column mapping
    *na_values*.
    """
    return PD.read_fwf(omni_fname,
                       colspecs=colspecs,
                       header=None,
                       names=names,
                       na_values=na_values)
