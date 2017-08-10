import decimal

"""
SI unit prefix names, from yocto to yotta.
"""
PREFIX = u'yzafpn\u03BCm kMGTPEZY'

"""
Power 10 exponential corresponding to the min/max in *PREFIX*.
"""
SHIFT = decimal.Decimal('1E24')

def sistr(x, baseunit=''):
    """
    Return the string representation of *x* using the SI convention,
    appending *baseunit* if it is provided.

    Source: https://stackoverflow.com/questions/29627796/pretty-printing-physical-quantities-with-automatic-scaling-of-si-prefixes
    """
    if x == 0:
        return '0 ' + baseunit
    d = (decimal.Decimal(str(x)) * SHIFT).normalize()
    m, e = d.to_eng_string().split('E')
    return (m + ' ' + PREFIX[int(e) // 3] + baseunit).encode('utf-8')


if __name__ == '__main__':
    print(sistr(10100, 'g'))
    print(sistr(10100e-12, 'g'))
    print(sistr(12))
    print(sistr(1200))
    print(sistr(0))
    print(sistr(1e-6))
