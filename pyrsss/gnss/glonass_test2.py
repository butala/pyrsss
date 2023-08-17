import scipy.constants as const
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

#from constants import M_TO_TECU, NS_TO_TECU
from .rinex_new import RinexDump
from ..ionex.read_ionex import read_header
from .constants import K, F_GLO_1, F_GLO_2, F_GLO_1_DELTA, F_GLO_2_DELTA
from .glonass import GLONASS_Status


if __name__ == '__main__':
    dump_fname = '/tmp/algo0010.14o.glonass.dump'

    rinex_dump = RinexDump.load(dump_fname, p1c1=False)

    ionex_code_fname = '/Users/butala/Documents/IARPA/absolute_tec/codg0010.14i'

    with open(ionex_code_fname) as fid:
        _, sat_biases_code, stn_biases_code = read_header(fid)

    glonass_status = GLONASS_Status()

    fig = plt.figure(figsize=(11, 8.5))

    for sat in sorted(set(rinex_dump.sat.values)):
        info = glonass_status(int(sat[1:]), rinex_dump.gps_time[0])
        k = info.freq
        f1 = F_GLO_1 + k * F_GLO_1_DELTA
        f2 = F_GLO_2 + k * F_GLO_2_DELTA
        # f1 = F_GLO_1
        # f2 = F_GLO_2
        TECU_TO_NS =  K * 10**25 / const.c * (f1**2 - f2**2) / (f1**2 * f2**2)
        # TECU_TO_NS2 =  K * 10**25 / const.c * (F_GLO_1**2 - F_GLO_2**2) / (F_GLO_1**2 * F_GLO_2**2)
        NS_TO_TECU = 1 / TECU_TO_NS
        TECU_TO_KM = TECU_TO_NS * (const.c / (1e9 * 1e3))
        TECU_TO_M = TECU_TO_KM * 1e3
        M_TO_TECU = 1 / TECU_TO_M
        # stn_bias = stn_biases_code['GLONASS']['ALGO'][0] / TECU_TO_NS2
        stn_bias = stn_biases_code['GLONASS']['ALGO'][0] * NS_TO_TECU
        dump_i = rinex_dump[rinex_dump.sat == sat]
        dt = dump_i.gps_time.values
        sat_bias = sat_biases_code['GLONASS'][int(sat[1:])][0] * NS_TO_TECU
        # sat_bias = sat_biases_code['GLONASS'][int(sat[1:])][0] / TECU_TO_NS2
        # sat_bias = stn_bias = 0
        P_I = (dump_i.P2.values - dump_i.P1.values) * M_TO_TECU + sat_bias + stn_bias

        # stn_bias_ns = stn_biases_code['GLONASS']['ALGO'][0]
        # sat_bias_ns = sat_biases_code['GLONASS'][int(sat[1:])][0]

        # stn_bias_s = stn_bias_ns * 1e-9
        # sat_bias_s = sat_bias_ns * 1e-9

        # stn_bias_m = stn_bias_s * const.c
        # sat_bias_m = sat_bias_s * const.c

        # P_I_m = (dump_i.P2.values - dump_i.P1.values) + sat_bias_m + stn_bias_m

        # fig.clf()
        # plt.plot_date(dt, P_I, marker='x', color='r')
        plt.plot_date(dt, P_I, marker='.', color='r', ms=1)
        # plt.plot_date(dt, P_I_m, marker='.', color='r', ms=1)

    # plt.ylim(-5, 100)
    plt.savefig('/tmp/glonass_test2.pdf'.format(sat), bbox_inches='tight')
