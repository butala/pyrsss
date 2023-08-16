import scipy as sp


def minute_interval_filter_firwin():
    """
    Synthesize and return impulse response of minute interval filter
    for analysis of magnetometer data. Design a filter using the
    window method. The filter is designed to stop from 0 to 0.75 mHz
    and pass from 1 mHz to Nyquist.

    For reference: this is the filter used to process the data for the
    Space Weather 2017 paper.

    fs_1m = 1./60
    bands_1m = [0, 1e-3 - 0.25e-3, 1e-3, fs_1m/2]
    desired_1m = [0, 1]

    _ = fir_response(h_1m, bands_1m, desired_1m, Hz=fs_1m)

    FIR filter is linear phase type 1

    band 0 abs deviations from 0:
    min = 4.466e-05  (db=-87.000846)
    med = -4.218e-03
    max = 3.027e-02 (db=-30.379491

    band 1 abs deviations from 1:
    min = 4.412e-06  (db=-107.106525)
    med = 2.398e-03
    std = 4.059e-03
    max = 2.982e-02 (db=-30.509428
    """
    Ts = 60
    fs = 1. / Ts
    fn = fs / 2
    width = 0.25e-3
    cutoff = 1e-3
    ripple = 30
    numtaps, beta = sp.signal.kaiserord(ripple, width / fn)
    if numtaps % 2 == 0:
        numtaps += 1
    return sp.signal.firwin(numtaps,
                            cutoff - width/2,
                            width=width,
                            pass_zero=False,
                            fs=fs)


def minute_interval_filter(N_remez=201):
    """
    Synthesize and return impulse response of minute interval filter
    for analysis of magnetometer data. Design a length *N_remez*
    min-max optimal filter. The filter is designed to stop from 0 to
    0.75 mHz and pass from 1 mHz to Nyquist.

    fs_1m = 1./60
    bands_1m = [0, 1e-3 - 0.25e-3, 1e-3, fs_1m/2]
    desired_1m = [0, 1]

    _ = fir_response(h_1m, bands_1m, desired_1m, Hz=fs_1m)

    FIR filter is linear phase type 1

    band 0 abs deviations from 0:
    min = 1.487e-05  (db=-96.553225)
    med = -1.102e-03
    max = 1.569e-03 (db=-56.089238

    band 1 abs deviations from 1:
    min = 1.018e-06  (db=-119.848775)
    med = -1.066e-05
    std = 4.820e-04
    max = 1.571e-03 (db=-56.074809
    """
    return sp.signal.remez(N_remez,
                           [0, 0.75e-3, 1e-3, 1./60/2],
                           [0, 1],
                           [1, 1],
                           Hz=1./60)


def second_interval_filter(N_remez=4001):
    """
    Synthesize and return impulse response of second interval filter
    for analysis of magnetometer data. Design a length *N_remez*
    min-max optimal filter. The filter is designed to stop from 0 to
    0.25 mHz, pass from 1 to 100 mHz, and stop from 101 mHz to
    Nyquist.


    fs_1m = 1.
    bands_1m = [0, 1e-3 - 0.25e-3,
                1e-3, 100e-3,
                101e-3, fs_1m/2]
    desired_1m = [0, 1, 0]

    _ = fir_response(h_1m, bands_1m, desired_1m, Hz=fs_1m)

    FIR filter is linear phase type 1

    band 0: 0 Hz -- 750 \muHz
    abs deviations from 0 statistics:
    min = 6.872e-05  (db=-83.257739)
    med = -4.906e-02
    max = 7.486e-01 (db=-2.514808)

    band 1: 1 mHz -- 100 mHz
    abs deviations from 1 statistics:
    min = 2.323e-07  (db=-132.678419)
    med = -6.093e-07
    std = 5.871e-04
    max = 2.349e-03 (db=-52.584106)

    band 2: 101 mHz -- 500 mHz
    abs deviations from 0 statistics:
    min = 1.430e-07  (db=-136.893586)
    med = -1.353e-03
    max = 2.853e-03 (db=-50.895031)
    """
    return sp.signal.remez(N_remez,
                           [0, 0.25e-3, 1e-3, 100e-3, 101e-3, 0.5],
                           [0, 1, 0],
                           [1, 1, 1],
                           Hz=1)


if __name__ == '__main__':
    from pyrsss.signal.lfilter import fir_response

    h_1m = minute_interval_filter_firwin()

    fs_1m = 1.
    bands_1m = [0, 1e-3 - 0.25e-3,
                1e-3, 100e-3,
                101e-3, fs_1m/2]
    desired_1m = [0, 1, 0]

    _ = fir_response(h_1m, bands_1m, desired_1m, Hz=fs_1m)
