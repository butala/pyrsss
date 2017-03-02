from collections import OrderedDict, namedtuple

"""
Information found at http://vmo.igpp.ucla.edu, looking for MEASURE
in the "Find Data" section, and then looking at the XML record
associated with each observatory.
"""


class Info(namedtuple('Info', 'lat lon')):
    pass


MEASURE_INFO = OrderedDict([('APL', Info(39.170, -76.880)),
                            ('FIT', Info(28.070, -80.950)),
                            ('CLK', Info(44.700, -75.000)),
                            ('DSO', Info(36.250, -81.400)),
                            ('GTF', Info(43.620, -71.950)),
                            ('JAX', Info(30.350, -81.600)),
                            ('MSH', Info(42.600, -71.480)),
                            ('USC', Info(33.340, -81.460))])
